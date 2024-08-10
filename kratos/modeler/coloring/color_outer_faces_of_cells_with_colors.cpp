//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "color_outer_faces_of_cells_with_colors.h"

namespace Kratos {

ColorOuterFacesOfCellsWithColors::ColorOuterFacesOfCellsWithColors(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    VoxelMesherColoring(rModeler, ColoringParameters)
{}


void ColorOuterFacesOfCellsWithColors::Apply() const
{
    Parameters parameters = GetParameters();
    CartesianMeshColors& r_colors = GetMeshColors();

    const int interface_color = parameters["color"].GetInt();
    const int cell_color = parameters["cell_color"].GetInt();

    array_1d<std::size_t, 3> number_of_cells;
    for(int i = 0 ; i < 3 ; i++){
        number_of_cells[i] = GetKeyPlanes(i).size() - 1;
    }

    array_1d<std::size_t, 3> cell_indices;
    for (cell_indices[2] = 0; cell_indices[2] < number_of_cells[2]; cell_indices[2]++) {
        for (cell_indices[1] = 0; cell_indices[1] < number_of_cells[1]; cell_indices[1]++) {
            for (cell_indices[0] = 0; cell_indices[0] < number_of_cells[0]; cell_indices[0]++) {
                if(std::lround(r_colors.GetElementalColor(cell_indices[0],cell_indices[1],cell_indices[2])) == cell_color){
                    ApplyColorIfOuterFace(interface_color, cell_color, cell_indices);
                }
            }
        }
    }
}


void ColorOuterFacesOfCellsWithColors::ApplyColorIfOuterFace(
    const int InterfaceColor,
    const int CellColor,
    const array_1d<std::size_t, 3>& rCellIndices) const
{
    CartesianMeshColors& r_colors = GetMeshColors();
    auto& r_faces_color = r_colors.GetElementalFaceColor(rCellIndices[0],rCellIndices[1],rCellIndices[2]);
    for(int i_direction = 0 ; i_direction < 3 ; i_direction++){
        const auto& r_key_planes = GetKeyPlanes(i_direction);
        if(rCellIndices[i_direction] == 0) {  // It is the first cell and has outer face
            r_faces_color[i_direction] = InterfaceColor;
        } else {
            array_1d<std::size_t, 3> neighbor_cell_indices = rCellIndices;
            neighbor_cell_indices[i_direction]--;
            const int neigh_cell_color = std::lround(r_colors.GetElementalColor(neighbor_cell_indices[0],neighbor_cell_indices[1],neighbor_cell_indices[2]));
            if(neigh_cell_color != CellColor){
                r_faces_color[i_direction] = InterfaceColor;
            }
        }
        if(rCellIndices[i_direction] == (r_key_planes.size() - 2)) {  // It is the last cell has outer face
            r_faces_color[i_direction + 3] = InterfaceColor;
        } else {
            array_1d<std::size_t, 3> neighbor_cell_indices = rCellIndices;
            neighbor_cell_indices[i_direction]++;
            const int neigh_cell_color = std::lround(r_colors.GetElementalColor(neighbor_cell_indices[0],neighbor_cell_indices[1],neighbor_cell_indices[2]));
            if(neigh_cell_color != CellColor){
                r_faces_color[i_direction + 3] = InterfaceColor;
            }
        }
    }
}


}