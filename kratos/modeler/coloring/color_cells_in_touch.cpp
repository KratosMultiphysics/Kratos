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
#include "color_cells_in_touch.h"

namespace Kratos {

ColorCellsInTouch::ColorCellsInTouch(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    VoxelMesherColoring(rModeler, ColoringParameters)
{}


void ColorCellsInTouch::Apply() const
{
    auto parameters = GetParameters();
    auto& r_model_part = GetModelPart(parameters["model_part_name"].GetString());

    const int inside_color = parameters["color"].GetInt();
    const int outside_color = parameters["outside_color"].GetInt();
    const std::string input_entities = parameters["input_entities"].GetString();

    if (input_entities == "elements") {
        for(auto& r_geometrical_object : r_model_part.Elements()) {
            const auto& r_geometry = r_geometrical_object.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToCellsInTouchWithGeometry(r_geometry, inside_color, outside_color);
        }
    } else if(input_entities == "conditions") {
        for(auto& r_geometrical_object : r_model_part.Conditions()) {
            const auto& r_geometry = r_geometrical_object.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToCellsInTouchWithGeometry(r_geometry, inside_color, outside_color);
        }
    } else if(input_entities == "geometries") {
        for(const auto& r_geometry : r_model_part.Geometries()) {
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToCellsInTouchWithGeometry(r_geometry, inside_color, outside_color);
        }
    } else {
        KRATOS_ERROR << "The input_entities  " << parameters["input_entities"] << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
    }
}


void ColorCellsInTouch::ApplyColorToCellsInTouchWithGeometry(
    const Geometry<Node>& rGeometry,
    int InsideColor,
    int OutsideColor) const
{
    auto& r_colors = GetMeshColors();

    array_1d<std::size_t, 3> min_position(3,0);
    array_1d<std::size_t, 3> max_position(3,0);
    r_colors.CalculateOuterMinMaxNodePositions(rGeometry, min_position, max_position);
    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
        for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
            for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                Point cell_min_point = r_colors.GetPoint(i,j,k);
                Point cell_max_point = r_colors.GetPoint(i+1,j+1,k+1);
                if(rGeometry.HasIntersection(cell_min_point,cell_max_point)){
                    r_colors.GetElementalColor(i,j,k) = InsideColor;
                }
            }
        }
    }
}

}