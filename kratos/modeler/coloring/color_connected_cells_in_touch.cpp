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
#include "color_connected_cells_in_touch.h"

namespace Kratos {

ColorConnectedCellsInTouch::ColorConnectedCellsInTouch(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    VoxelMesherColoring(rModeler, ColoringParameters)
{}


void ColorConnectedCellsInTouch::Apply() const
{
    Parameters parameters = GetParameters();
    ModelPart& r_skin_model_part = GetModelPart(parameters["model_part_name"].GetString());

    const int inside_color = parameters["color"].GetInt();
    const int cell_color = parameters["cell_color"].GetInt();
    const std::string input_entities = parameters["input_entities"].GetString();

    if(input_entities == "elements") {
        for(auto& r_geometrical_object : r_skin_model_part.Elements()) {
            auto& r_geometry = r_geometrical_object.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToConnectedCellsInTouchWithGeometry(r_geometry, inside_color, cell_color);
        }
    } else if(input_entities == "conditions") {
        for(auto& r_geometrical_object : r_skin_model_part.Conditions()) {
            auto& r_geometry = r_geometrical_object.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToConnectedCellsInTouchWithGeometry(r_geometry, inside_color, cell_color);
        }
    } else if(input_entities == "geometries") {
        for(auto& r_geometry : r_skin_model_part.Geometries()) {
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToConnectedCellsInTouchWithGeometry(r_geometry, inside_color, cell_color);
        }
    } else{
        KRATOS_ERROR << "The input_entities " << input_entities << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
    }
}


void ColorConnectedCellsInTouch::ApplyColorToConnectedCellsInTouchWithGeometry(
    const Geometry<Node>& rGeometry,
    const int InsideColor,
    const int CellColor
    ) const
{
    array_1d<std::size_t, 3> min_position(3,0);
    array_1d<std::size_t, 3> max_position(3,0);
    CartesianMeshColors& r_colors = GetMeshColors();
    r_colors.CalculateOuterMinMaxNodePositions(rGeometry, min_position, max_position);
    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
        for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
            for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                const int cell_color = std::lround(r_colors.GetElementalColor(i,j,k));
                if (cell_color == CellColor) {
                    Point cell_min_point = r_colors.GetPoint(i,j,k);
                    Point cell_max_point = r_colors.GetPoint(i+1,j+1,k+1);
                    if(rGeometry.HasIntersection(cell_min_point,cell_max_point)){
                        ColorConnectedCellsToThisCell(i,j,k,InsideColor, CellColor, r_colors);
                    }
                }
            }
        }
    }
}


void ColorConnectedCellsInTouch::ColorConnectedCellsToThisCell(
    const std::size_t I,
    const std::size_t J,
    const std::size_t K,
    const int InsideColor,
    const int CellColor,
    CartesianMeshColors& rColors) const
{
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> stack;
    stack.emplace_back(I,J,K);
    const std::size_t i_size = rColors.GetNodalCoordinates(0).size();
    const std::size_t j_size = rColors.GetNodalCoordinates(1).size();
    const std::size_t k_size = rColors.GetNodalCoordinates(2).size();
    while (!stack.empty()) {
        auto& indices = stack.back();
        const auto i = std::get<0>(indices);
        const auto j = std::get<1>(indices);
        const auto k = std::get<2>(indices);
        stack.pop_back();
        rColors.GetElementalColor(i,j,k) = InsideColor;
        if(i > 0){
            if(std::lround(rColors.GetElementalColor(i-1, j, k)) == CellColor){
                stack.emplace_back(i-1, j, k);
            }
        }
        if(j > 0){
            if(std::lround(rColors.GetElementalColor(i, j-1, k)) == CellColor){
                stack.emplace_back(i, j-1, k);
            }
        }
        if(k > 0){
            if(std::lround(rColors.GetElementalColor(i, j, k-1)) == CellColor){
                stack.emplace_back(i, j, k-1);
            }
        }
        if(i+1 < i_size){
            if(std::lround(rColors.GetElementalColor(i+1, j, k)) == CellColor){
                stack.emplace_back(i+1, j, k);
            }
        }
        if(j+1 < j_size){
            if(std::lround(rColors.GetElementalColor(i, j+1, k)) == CellColor){
                stack.emplace_back(i, j+1, k);
            }
        }
        if(k+1 < k_size){
            if(std::lround(rColors.GetElementalColor(i, j, k+1)) == CellColor){
                stack.emplace_back(i, j, k+1);
            }
        }
    }
}

}
