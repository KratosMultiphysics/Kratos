//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher

#if !defined(KRATOS_MAPPING_UTILITIES_H_INCLUDED )
#define  KRATOS_MAPPING_UTILITIES_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos
{
double ComputeSearchRadius(ModelPart& model_part) {
    double search_safety_factor = 1.2;
    double max_element_size = 0.0;

    // Loop through each edge of a geometrical entity ONCE
    for (auto& condition : model_part.GetCommunicator().LocalMesh().Conditions()) {
        for (std::size_t i = 0; i < (condition.GetGeometry().size() - 1); ++i) {
            double node_1_x = condition.GetGeometry()[i].X();
            double node_1_y = condition.GetGeometry()[i].Y();
            double node_1_z = condition.GetGeometry()[i].Z();

            for (std::size_t j = i + 1; j < condition.GetGeometry().size(); ++j) {
                double node_2_x = condition.GetGeometry()[j].X();
                double node_2_y = condition.GetGeometry()[j].Y();
                double node_2_z = condition.GetGeometry()[j].Z();

                double edge_length = sqrt(pow(node_1_x - node_2_x , 2) + pow(node_1_y - node_2_y , 2) + pow(node_1_z - node_2_z , 2));
                max_element_size = std::max(max_element_size, edge_length);
            }
        }
    }
    model_part.GetCommunicator().MaxAll(max_element_size); // Compute the maximum among the partitions
    return max_element_size * search_safety_factor;
}

}  // namespace Kratos.

#endif // KRATOS_MAPPING_UTILITIES_H_INCLUDED  defined
