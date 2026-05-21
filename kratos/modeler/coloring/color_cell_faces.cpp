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
#include "color_cell_faces.h"

namespace Kratos {

ColorCellFaces::ColorCellFaces(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    VoxelMesherColoring(rModeler, ColoringParameters)
{}


void ColorCellFaces::Apply() const
{
    auto parameters = GetParameters();
    auto& r_model_part = GetModelPart(parameters["model_part_name"].GetString());
    CartesianMeshColors& r_colors = GetMeshColors();

    const int outside_color = parameters["outside_color"].GetInt();
    const int interface_color = parameters["color"].GetInt();
    const int cell_color = parameters["cell_color"].GetInt();
    const std::string input_entities = parameters["input_entities"].GetString();

    array_1d< std::size_t, 3 > min_ray_position{0, 0, 0};
    array_1d< std::size_t, 3 > max_ray_position{0, 0, 0};

    r_colors.CalculateMinMaxCenterOfElementPositions(r_model_part.Nodes(), min_ray_position, max_ray_position);
    r_colors.InitializeRays(min_ray_position, max_ray_position, "center_of_elements");

    if(input_entities == "elements") {
        for(auto& element : r_model_part.Elements()) {
            Element::GeometryType& r_geometry = element.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            r_colors.AddGeometry(r_geometry, false);
        }
    } else if(input_entities == "conditions") {
        for(auto& condition : r_model_part.Conditions()) {
            auto& r_geometry = condition.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            r_colors.AddGeometry(r_geometry, false);
        }
    } else if(input_entities == "geometries") {
        for(auto& r_geometry : r_model_part.Geometries()) {
            CheckGeometryType(r_geometry.GetGeometryType());
            r_colors.AddGeometry(r_geometry, false);
        }
    }

    r_colors.CalculateElementalFaceColors(min_ray_position, max_ray_position, interface_color, outside_color, cell_color);
}

}