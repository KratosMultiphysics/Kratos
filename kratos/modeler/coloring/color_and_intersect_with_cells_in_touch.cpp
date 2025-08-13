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
#include "color_and_intersect_with_cells_in_touch.h"

namespace Kratos {

ColorAndIntersectWithCellsInTouch::ColorAndIntersectWithCellsInTouch(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    VoxelMesherColoring(rModeler, ColoringParameters)
{}


void ColorAndIntersectWithCellsInTouch::Apply() const
{
    auto parameters = GetParameters();
    auto& r_model_part = GetModelPart(parameters["model_part_name"].GetString());

    const int inside_color = parameters["color"].GetInt();
    const std::string input_entities = parameters["input_entities"].GetString();

    auto& r_skin_intersections = GetIntersections(inside_color);

    if (input_entities == "elements") {
        for(auto& r_geometrical_object : r_model_part.Elements()) {
            auto& r_geometry = r_geometrical_object.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToCellsInTouchWithGeometryAndComputeIntersection(r_geometry, inside_color, r_skin_intersections);
        }
    } else if(input_entities == "conditions") {
        for(auto& r_geometrical_object : r_model_part.Conditions()) {
            auto& r_geometry = r_geometrical_object.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToCellsInTouchWithGeometryAndComputeIntersection(r_geometry, inside_color, r_skin_intersections);
        }
    } else if(input_entities == "geometries") {
        for(auto& r_geometry : r_model_part.Geometries()) {
            CheckGeometryType(r_geometry.GetGeometryType());
            ApplyColorToCellsInTouchWithGeometryAndComputeIntersection(r_geometry, inside_color, r_skin_intersections);
        }
    } else {
        KRATOS_ERROR << "The input_entities  " << parameters["input_entities"] << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
    }
}


Parameters ColorAndIntersectWithCellsInTouch::GetDefaultParameters() const {
    return Parameters(R"({
        "type" : "Undefined_coloring_type",
        "model_part_name": "Undefined",
        "color": -1,
        "cell_color": -1,
        "input_entities": "elements",
        "outside_color": 0
    })");
}


void ColorAndIntersectWithCellsInTouch::ApplyColorToCellsInTouchWithGeometryAndComputeIntersection(
    const Element::GeometryType& rGeometry,
    int InsideColor,
    Internals::SkinIntersection& rSkinIntersection) const
{
    auto& r_colors = GetMeshColors();

    auto p_work_geometry = rSkinIntersection.ConvertToWorkGeometry(rGeometry);
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

                    rSkinIntersection.ComputeCutInsideCell(p_work_geometry, cell_min_point, cell_max_point, i, j, k);
                }
            }
        }
    }
}

}