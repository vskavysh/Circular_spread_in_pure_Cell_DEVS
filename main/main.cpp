#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <cadmium/modeling/celldevs/grid/coupled.hpp>
#include <cadmium/simulation/logger/csv.hpp>
#include <cadmium/simulation/root_coordinator.hpp>

#include "addGridCell.hpp"

using namespace cadmium;
using namespace cadmium::celldevs;

namespace fs = std::filesystem;

//This derives the name of the log file from the name of the config file
static std::string derive_log_filename(const std::string& config_path_str) {
    fs::path config_path(config_path_str);

    std::string stem = config_path.stem().string();
    const std::string config_suffix = "_config";

    if (stem.size() >= config_suffix.size() &&
        stem.compare(stem.size() - config_suffix.size(), config_suffix.size(), config_suffix) == 0) {
        stem.replace(stem.size() - config_suffix.size(), config_suffix.size(), "_log");
    } else {
        stem += "_log";
    }

    fs::path out = config_path.parent_path() / (stem + ".csv");
    return out.string();
}

//This removes the "sep=;" line at the beginning of the generated log file. This line causes an error in DEVS Simulator Viewer
static void remove_sep_line_if_present(const std::string& file_path) {
    std::ifstream in(file_path);
    if (!in.is_open()) {
        std::cerr << "Warning: could not reopen log file: " << file_path << "\n";
        return;
    }

    std::string first_line;
    if (!std::getline(in, first_line)) {
        return;
    }

    if (first_line != "sep=;") {
        return;
    }

    std::string remainder;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (!first) remainder += "\n";
        remainder += line;
        first = false;
    }
    in.close();

    std::ofstream out(file_path, std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "Warning: could not rewrite log file: " << file_path << "\n";
        return;
    }
    out << remainder;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0]
                  << " SCENARIO_CONFIG.json [MAX_SIMULATION_TIME_MS default: 20.0]"
                  << " [LOG_FILE default: derived from config]\n";
        return -1;
    }

    const std::string config_file = argv[1];
    const double sim_time = (argc > 2) ? std::stod(argv[2]) : 20.0;
    const std::string log_file = (argc > 3) ? argv[3] : derive_log_filename(config_file);

    try {
        auto model = std::make_shared<GridCellDEVSCoupled<brainWaveState, double>>(
            "brain_wave_grid",
            addGridCell,
            config_file
        );

        model->buildModel();

        auto rootCoordinator = RootCoordinator(model);
        rootCoordinator.setLogger<CSVLogger>(log_file, ";");

        rootCoordinator.start();
        rootCoordinator.simulate(sim_time);
        rootCoordinator.stop();

        remove_sep_line_if_present(log_file);

        std::cout << "Simulation finished.\n";
        std::cout << "Config: " << config_file << "\n";
        std::cout << "Log:    " << log_file << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Simulation failed: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
