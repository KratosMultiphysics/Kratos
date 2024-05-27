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
#include "includes/dem_variables.h"
#include "color_cells_with_center_in_sphere_around_nodes.h"

namespace Kratos {

ColorCellsWithCenterInSphereArounNodes::ColorCellsWithCenterInSphereArounNodes(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    VoxelMesherColoring(rModeler, ColoringParameters)
{}


void ColorCellsWithCenterInSphereArounNodes::Apply() const
{
    const ModelPart& r_input_model_part = GetInputModelPart();
    CartesianMeshColors& r_colors = GetMeshColors();
    Parameters parameters = GetParameters();
    const int color = parameters["color"].GetInt();

    for(const auto& i_node : r_input_model_part.Nodes()){
        array_1d< std::size_t, 3 > min_position;
        array_1d< std::size_t, 3 > max_position;
        double radius = i_node.GetSolutionStepValue(RADIUS);
        double radius2 = radius * radius; // storing radius**2 to avoid sqrt in comparison

        for(int i = 0; i < 3; i++ ) {
            min_position[i] = r_colors.CalculateCenterOfElementPosition(i_node[i] - radius, i);
            max_position[i] = r_colors.CalculateCenterOfElementPosition(i_node[i] + radius, i) + 1;
        }

        for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
            for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    Point cell_center = r_colors.GetCenterOfElement(i,j,k);
                    array_1d<double, 3> distance_vector = i_node.Coordinates() -  cell_center.Coordinates();
                    if(inner_prod(distance_vector, distance_vector) < radius2){
                        r_colors.GetElementalColor(i,j,k) = color;
                    }
                }
            }
        }
    }
}

}