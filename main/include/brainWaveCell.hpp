#ifndef BRAIN_WAVE_CELL_HPP
#define BRAIN_WAVE_CELL_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <unordered_map>
#include <vector>

#include <cadmium/modeling/celldevs/grid/cell.hpp>
#include <cadmium/modeling/celldevs/grid/config.hpp>

#include "brainWaveState.hpp"

using namespace cadmium::celldevs;

class brainWave : public GridCell<brainWaveState, double> {
public:
    explicit brainWave(
        const std::vector<int>& id,
        const std::shared_ptr<const GridCellConfig<brainWaveState, double>>& config
    )
        : GridCell<brainWaveState, double>(id, config),
          myId(id),
          gen(std::random_device{}()),
          uni01(0.0, 1.0) {}

    [[nodiscard]] brainWaveState localComputation(
        brainWaveState state,
        const std::unordered_map<std::vector<int>, NeighborData<brainWaveState, double>>& neighborhood
    ) const override {
        brainWaveState next = state;
        next.current_outgoing.clear();

        if (next.is_neuron == 1) {
            advanceNeuron(next);
            return next;
        }

        auto accepted = gatherAcceptedArrivals(next, myId, neighborhood);

        if (accepted.empty()) {
            next.potential = REST_MV;
            next.sigma = std::numeric_limits<double>::infinity();
            return next;
        }

        // For each neuron id, keep exactly one message with the highest sequence.
        accepted = keepOnlyLatestOnePerNeuron(accepted);

        if (accepted.empty()) {
            next.potential = REST_MV;
            next.sigma = std::numeric_limits<double>::infinity();
            return next;
        }

        next.potential = computeAveragePotential(accepted);
        next.current_outgoing = expandAll(accepted);

        if (next.current_outgoing.empty()) {
            next.sigma = std::numeric_limits<double>::infinity();
        } else {
            next.sigma = averageDelay(next.current_outgoing);
        }

        return next;
    }

    [[nodiscard]] double outputDelay(const brainWaveState& state) const override {
        return state.sigma;
    }

private:
    static constexpr double REST_MV = -70.0;
    static constexpr double SPIKE_PEAK_MV = 30.0;
    static constexpr double DELTA_TAU_MS = 0.01;
    static constexpr double REFRACTORY_MS = 1.9;
    static constexpr double TRIANGLE_MS = 2.0;
    static constexpr double DISTANCE_SLOWDOWN = 10.0;
    static constexpr double EPS = 1e-12;

    //// Round all delays to nearest DELTA_TAU_MS / 10
    //static constexpr double DELAY_QUANTUM_MS = DELTA_TAU_MS / 10.0;
    // Round all delays to nearest DELTA_TAU_MS / 1
    static constexpr double DELAY_QUANTUM_MS = DELTA_TAU_MS / 1.0;

    std::vector<int> myId;
    mutable std::mt19937 gen;
    mutable std::uniform_real_distribution<double> uni01;

    static double quantizeDelay(double dt) {
        if (!std::isfinite(dt)) return dt;
        if (dt <= 0.0) return 0.0;
        return DELAY_QUANTUM_MS * std::round(dt / DELAY_QUANTUM_MS);
    }

    static double distFromOrigin(int m, int n) {
        return std::sqrt(static_cast<double>(m * m + n * n));
    }

    static double triangularVoltage(double age_ms) {
        if (age_ms < 0.0 || age_ms > TRIANGLE_MS) return REST_MV;
        if (age_ms <= 1.0) return REST_MV + 100.0 * age_ms;
        return SPIKE_PEAK_MV - 100.0 * (age_ms - 1.0);
    }

    static double messageDelay(const WaveMsg& msg) {
        const double raw =
            DELTA_TAU_MS *
            (distFromOrigin(msg.m_intended, msg.n_intended) -
             distFromOrigin(msg.m, msg.n));

        return quantizeDelay(raw);
    }

    static bool sameMessage(const WaveMsg& a, const WaveMsg& b) {
        return std::abs(a.V0 - b.V0) <= EPS &&
               a.m == b.m &&
               a.n == b.n &&
               a.m_intended == b.m_intended &&
               a.n_intended == b.n_intended &&
               a.n_sequence == b.n_sequence &&
               a.n_neuron == b.n_neuron;
    }

    static void addUnique(std::vector<WaveMsg>& out, const WaveMsg& msg) {
        for (const auto& kept : out) {
            if (sameMessage(kept, msg)) return;
        }
        out.push_back(msg);
    }

    void advanceNeuron(brainWaveState& state) const {
        state.tick += 1;
        state.current_outgoing.clear();

        if (state.firing == 1) {
            state.spike_age_ms += DELTA_TAU_MS;

            if (state.spike_age_ms > TRIANGLE_MS) {
                state.firing = 0;
                state.spike_age_ms = 0.0;
                state.potential = REST_MV;
                state.refractory_remaining_ms = REFRACTORY_MS;
            } else {
                state.potential = triangularVoltage(state.spike_age_ms);
            }
        } else {
            state.potential = REST_MV;

            if (state.refractory_remaining_ms > 0.0) {
                state.refractory_remaining_ms =
                    std::max(0.0, state.refractory_remaining_ms - DELTA_TAU_MS);
            } else {
                const bool fire_now = (uni01(gen) < state.p_fire);
                if (fire_now) {
                    state.firing = 1;
                    state.spike_age_ms = 0.0;
                    state.potential = triangularVoltage(0.0);
                }
            }
        }

        if (state.firing == 1) {
            const double Vemit = triangularVoltage(state.spike_age_ms);

            if (Vemit > REST_MV + EPS) {
                const std::uint64_t n_sequence = state.tick;

                state.current_outgoing = {
                    WaveMsg(Vemit, 0, 0,  1,  0, n_sequence, state.n_neuron),
                    WaveMsg(Vemit, 0, 0, -1,  0, n_sequence, state.n_neuron),
                    WaveMsg(Vemit, 0, 0,  0,  1, n_sequence, state.n_neuron),
                    WaveMsg(Vemit, 0, 0,  0, -1, n_sequence, state.n_neuron)
                };
            }
        }

        state.sigma = DELTA_TAU_MS;
    }

    static bool alreadyProcessedByCurrentState(const brainWaveState& state, const WaveMsg& candidate) {
        for (const auto& msg : state.current_outgoing) {
            if (msg.n_neuron == candidate.n_neuron &&
                msg.n_sequence >= candidate.n_sequence) {
                return true;
            }
        }
        return false;
    }

    std::vector<WaveMsg> gatherAcceptedArrivals(
        const brainWaveState& currentState,
        const std::vector<int>& currentId,
        const std::unordered_map<std::vector<int>, NeighborData<brainWaveState, double>>& neighborhood
    ) const {
        std::vector<WaveMsg> accepted;

        for (const auto& [neighborId, neighborData] : neighborhood) {
            if (!neighborData.state) continue;
            const auto& neighborState = *(neighborData.state);

            const int dx = currentId[0] - neighborId[0];
            const int dy = currentId[1] - neighborId[1];

            for (const auto& msg : neighborState.current_outgoing) {
                const int m_current = msg.m + dx;
                const int n_current = msg.n + dy;

                if (m_current != msg.m_intended || n_current != msg.n_intended) {
                    continue;
                }

                WaveMsg adjusted(
                    msg.V0,
                    msg.m,
                    msg.n,
                    msg.m_intended,
                    msg.n_intended,
                    msg.n_sequence,
                    msg.n_neuron
                );

                if (alreadyProcessedByCurrentState(currentState, adjusted)) {
                    continue;
                }

                addUnique(accepted, adjusted);
            }
        }

        return accepted;
    }

    static std::vector<WaveMsg> keepOnlyLatestOnePerNeuron(const std::vector<WaveMsg>& msgs) {
        std::vector<WaveMsg> result;

        for (const auto& msg : msgs) {
            bool found = false;

            for (auto& kept : result) {
                if (kept.n_neuron == msg.n_neuron) {
                    found = true;

                    if (msg.n_sequence > kept.n_sequence) {
                        kept = msg;
                    }
                    break;
                }
            }

            if (!found) {
                result.push_back(msg);
            }
        }

        return result;
    }

    static double computeAveragePotential(const std::vector<WaveMsg>& msgs) {
        if (msgs.empty()) return REST_MV;

        double sum = 0.0;
        for (const auto& msg : msgs) {
            const double r = distFromOrigin(msg.m_intended, msg.n_intended);
            if (r > 0.0) {
                const double d_eff = 1.0 + (r - 1.0) / DISTANCE_SLOWDOWN;
                sum += (msg.V0 + 70.0) / d_eff;
            }
        }

        double V = REST_MV + sum / static_cast<double>(msgs.size());

        if (std::abs(V) < EPS) V = 0.0;
        if (std::abs(V - REST_MV) < EPS) V = REST_MV;
        return V;
    }

    static std::vector<WaveMsg> childrenOf(const WaveMsg& accepted) {
        std::vector<WaveMsg> out;

        const double V0 = accepted.V0;
        const int m = accepted.m_intended;
        const int n = accepted.n_intended;
        const std::uint64_t n_sequence = accepted.n_sequence;
        const int n_neuron = accepted.n_neuron;

        if (m > 0) {
            addUnique(out, WaveMsg(V0, m, n, m + 1, n, n_sequence, n_neuron));
        } else if (m < 0) {
            addUnique(out, WaveMsg(V0, m, n, m - 1, n, n_sequence, n_neuron));
        }

        if (n > 0) {
            addUnique(out, WaveMsg(V0, m, n, m, n + 1, n_sequence, n_neuron));
        } else if (n < 0) {
            addUnique(out, WaveMsg(V0, m, n, m, n - 1, n_sequence, n_neuron));
        }

        // Special case only for (+/-1, 0)
        if (m == 1 && n == 0) {
            addUnique(out, WaveMsg(V0, 1, 0, 2, 0, n_sequence, n_neuron));
            addUnique(out, WaveMsg(V0, 1, 0, 1, 1, n_sequence, n_neuron));
            addUnique(out, WaveMsg(V0, 1, 0, 1, -1, n_sequence, n_neuron));
        } else if (m == -1 && n == 0) {
            addUnique(out, WaveMsg(V0, -1, 0, -2, 0, n_sequence, n_neuron));
            addUnique(out, WaveMsg(V0, -1, 0, -1, 1, n_sequence, n_neuron));
            addUnique(out, WaveMsg(V0, -1, 0, -1, -1, n_sequence, n_neuron));
        }

        return out;
    }

    static std::vector<WaveMsg> expandAll(const std::vector<WaveMsg>& accepted) {
        std::vector<WaveMsg> out;
        for (const auto& msg : accepted) {
            const auto spawned = childrenOf(msg);
            for (const auto& s : spawned) {
                addUnique(out, s);
            }
        }
        return out;
    }

    static double averageDelay(const std::vector<WaveMsg>& msgs) {
        if (msgs.empty()) {
            return std::numeric_limits<double>::infinity();
        }

        double sum = 0.0;
        std::size_t count = 0;
        for (const auto& msg : msgs) {
            const double dt = messageDelay(msg);
            if (dt > EPS) {
                sum += dt;
                count += 1;
            }
        }

        if (count == 0) {
            return DELTA_TAU_MS;
        }

        return quantizeDelay(sum / static_cast<double>(count));
    }
};

#endif
