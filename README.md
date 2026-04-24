# Circular wave propagation in pure Cell-DEVS

This repository contains a Cell-DEVS model implemented in Cadmium v2 for simulating the circular spread of a wave.

### Signal Propagation
- Local propagation via neighboring cells
- Approximate 1/r attenuation
- Signal superposition supported

## Repository Structure

Source code:
- main/main.cpp
- main/include/brainWaveState.hpp
- main/include/brainWaveCell.hpp
- main/include/addGridCell.hpp

Scenarios:
Located in scenarios/ folder (JSON config files)

Logs:
Generated as *_log.csv in scenarios/

## Requirements

- C++17
- CMake
- Cadmium v2
- nlohmann/json

## Build

./build_sim.sh

## Run

- bin/0_scenario_21x21_offcenter_defaultfiring_config.json

## Visualization

Use DEVS Simulation Viewer:
https://devssim.carleton.ca/cell-devs-viewer/

- approximate_circular_propagation.webm
is a video file showing the approximate circular wave propagation achieved. 

## Report

- Report_BCI_Reading_Writing_CellDEVS.pdf the Appendix details the various attempts made to create a circular wave spread in pure Cell-DEVS.

## Author

- Vladimir Skavysh
- Carleton University
