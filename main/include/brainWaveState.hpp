#ifndef BRAIN_WAVE_STATE_HPP
#define BRAIN_WAVE_STATE_HPP

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

#if __has_include(<nlohmann/json.hpp>)
#include <nlohmann/json.hpp>
#elif __has_include(<json.hpp>)
#include <json.hpp>
#else
#error "nlohmann/json.hpp not found."
#endif

struct WaveMsg {
    double V0;       // source neuron voltage at emission
    int m;           // sender cell x-offset from source neuron
    int n;           // sender cell y-offset from source neuron
    int m_intended;  // intended receiver x-offset from source neuron
    int n_intended;  // intended receiver y-offset from source neuron
    std::uint64_t n_sequence; // monotonically increasing pulse sample index
    int n_neuron;    // fixed neuron id from JSON

    WaveMsg()
        : V0(-70.0),
          m(0),
          n(0),
          m_intended(0),
          n_intended(0),
          n_sequence(0),
          n_neuron(-1) {}

    WaveMsg(double V0_,
            int m_,
            int n_,
            int m_intended_,
            int n_intended_,
            std::uint64_t n_sequence_,
            int n_neuron_)
        : V0(V0_),
          m(m_),
          n(n_),
          m_intended(m_intended_),
          n_intended(n_intended_),
          n_sequence(n_sequence_),
          n_neuron(n_neuron_) {}
};

inline bool operator!=(const WaveMsg& a, const WaveMsg& b) {
    return std::abs(a.V0 - b.V0) > 1e-12 ||
           a.m != b.m ||
           a.n != b.n ||
           a.m_intended != b.m_intended ||
           a.n_intended != b.n_intended ||
           a.n_sequence != b.n_sequence ||
           a.n_neuron != b.n_neuron;
}

struct brainWaveState {
    double potential;

    int firing;
    int is_neuron;
    double p_fire;
    int n_neuron;  // unique id for neurons, default -1 for passive tissue

    double spike_age_ms;
    double refractory_remaining_ms;

    double sigma;

    std::uint64_t tick;

    std::vector<WaveMsg> current_outgoing;

    brainWaveState()
        : potential(-70.0),
          firing(0),
          is_neuron(0),
          p_fire(0.02),
          n_neuron(-1),
          spike_age_ms(0.0),
          refractory_remaining_ms(0.0),
          sigma(std::numeric_limits<double>::infinity()),
          tick(0),
          current_outgoing() {}
};

inline bool operator!=(const brainWaveState& a, const brainWaveState& b) {
    if (std::abs(a.potential - b.potential) > 1e-12) return true;
    if (a.firing != b.firing) return true;
    if (a.is_neuron != b.is_neuron) return true;
    if (std::abs(a.p_fire - b.p_fire) > 1e-12) return true;
    if (a.n_neuron != b.n_neuron) return true;
    if (std::abs(a.spike_age_ms - b.spike_age_ms) > 1e-12) return true;
    if (std::abs(a.refractory_remaining_ms - b.refractory_remaining_ms) > 1e-12) return true;
    if (std::abs(a.sigma - b.sigma) > 1e-12) return true;
    if (a.tick != b.tick) return true;
    if (a.current_outgoing.size() != b.current_outgoing.size()) return true;

    for (std::size_t i = 0; i < a.current_outgoing.size(); ++i) {
        if (a.current_outgoing[i] != b.current_outgoing[i]) return true;
    }

    return false;
}

inline std::ostream& operator<<(std::ostream& os, const brainWaveState& s) {
    double out = s.potential;
    if (std::abs(out) < 1e-12) out = 0.0;
    if (std::abs(out + 70.0) < 1e-12) out = -70.0;
    os << out;
    return os;
}

inline void from_json(const nlohmann::json& j, brainWaveState& s) {
    s = brainWaveState();

    if (j.contains("potential")) j.at("potential").get_to(s.potential);
    if (j.contains("firing")) j.at("firing").get_to(s.firing);
    if (j.contains("is_neuron")) j.at("is_neuron").get_to(s.is_neuron);
    if (j.contains("p_fire")) j.at("p_fire").get_to(s.p_fire);
    if (j.contains("n_neuron")) j.at("n_neuron").get_to(s.n_neuron);
}

#endif
