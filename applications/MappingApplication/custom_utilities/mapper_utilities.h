//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_UTILITIES_H_INCLUDED)
#define  KRATOS_MAPPER_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
namespace MapperUtilities
{

/**
* @brief Assigning INTERFACE_EQUATION_IDs to the nodes, with and without MPI
* This function assigns the INTERFACE_EQUATION_IDs to the nodes, which
* act as EquationIds for the MappingMatrix. This work with and without MPI,
* in MPI a ScanSum is performed with the local number of nodes
* @param rModelPartCommunicator The Modelpart-Communicator to be used
* @author Philipp Bucher
*/
void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator);

inline int ComputeNumberOfNodes(ModelPart& rModelPart)
    {
        int num_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
        rModelPart.GetCommunicator().SumAll(num_nodes); // Compute the sum among the partitions
        return num_nodes;
    }

inline int ComputeNumberOfConditions(ModelPart& rModelPart)
{
    int num_conditions = rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
    rModelPart.GetCommunicator().SumAll(num_conditions); // Compute the sum among the partitions
    return num_conditions;
}

inline int ComputeNumberOfElements(ModelPart& rModelPart)
{
    int num_elements = rModelPart.GetCommunicator().LocalMesh().NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements); // Compute the sum among the partitions
    return num_elements;
}

inline double ComputeDistance(const array_1d<double, 3>& rCoords1,
                                const array_1d<double, 3>& rCoords2)
{
    return std::sqrt(std::pow(rCoords1[0] - rCoords2[0] , 2) +
                        std::pow(rCoords1[1] - rCoords2[1] , 2) +
                        std::pow(rCoords1[2] - rCoords2[2] , 2));
}

template <typename T>
inline double ComputeMaxEdgeLengthLocal(const T& rEntityContainer)
{
    double max_element_size = 0.0f;
    // Loop through each edge of a geometrical entity ONCE
    for (auto& r_entity : rEntityContainer)
    {
        for (std::size_t i = 0; i < (r_entity.GetGeometry().size() - 1); ++i)
        {
            for (std::size_t j = i + 1; j < r_entity.GetGeometry().size(); ++j)
            {
                double edge_length = ComputeDistance(r_entity.GetGeometry()[i].Coordinates(),
                                                        r_entity.GetGeometry()[j].Coordinates());
                max_element_size = std::max(max_element_size, edge_length);
            }
        }
    }
    return max_element_size;
}

inline double ComputeMaxEdgeLengthLocal(const ModelPart::NodesContainerType& rNodes)
{
    double max_element_size = 0.0f;
    // TODO modify loop such that it loop only once over the nodes
    for (auto& r_node_1 : rNodes)
    {
        for (auto& r_node_2 : rNodes)
        {
            double edge_length = ComputeDistance(r_node_1.Coordinates(),
                                                    r_node_2.Coordinates());
            max_element_size = std::max(max_element_size, edge_length);
        }
    }
    return max_element_size;
}

double ComputeSearchRadius(ModelPart& rModelPart, int EchoLevel);

inline double ComputeSearchRadius(ModelPart& rModelPart1, ModelPart& rModelPart2, int EchoLevel)
{
    double search_radius = std::max(ComputeSearchRadius(rModelPart1, EchoLevel),
                                    ComputeSearchRadius(rModelPart2, EchoLevel));
    return search_radius;
}

void CheckInterfaceModelParts(const int CommRank);

}  // namespace MapperUtilities.

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_H_INCLUDED  defined


