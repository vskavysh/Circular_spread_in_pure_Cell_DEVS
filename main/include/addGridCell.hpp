#ifndef ADD_GRID_CELL_HPP
#define ADD_GRID_CELL_HPP

#include <memory>
#include <stdexcept>

#include <cadmium/modeling/celldevs/grid/config.hpp>
#include <cadmium/modeling/celldevs/grid/cell.hpp>

#include "brainWaveCell.hpp"

using namespace cadmium::celldevs;

inline std::shared_ptr<GridCell<brainWaveState, double>> addGridCell(
    const coordinates& cellId,
    const std::shared_ptr<const GridCellConfig<brainWaveState, double>>& cellConfig) {
    const auto cellModel = cellConfig->cellModel;
    if (cellModel == "brainWave" || cellModel == "default") {
        return std::make_shared<brainWave>(cellId, cellConfig);
    }
    throw std::bad_typeid();
}

#endif
