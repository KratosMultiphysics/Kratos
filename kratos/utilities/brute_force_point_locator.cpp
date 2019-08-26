//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on work of Pablo Becker)
//

// System includes

// External includes

// Project includes
#include "brute_force_point_locator.h"

namespace Kratos
{

int BruteForcePointLocator::FindNode(const Point& rThePoint,
                                     const double DistanceThreshold) const
{
    int found_node_id = -1; // if no node is found this will be returned
    int local_node_found = 0;

    const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();
    const int global_num_nodes = r_data_comm.SumAll(static_cast<int>(mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes()));

    KRATOS_WARNING_IF("BruteForcePointLocator", global_num_nodes == 0) << r_data_comm << "No Nodes in ModelPart \"" << mrModelPart.Name() << std::endl;

    // note that this cannot be omp bcs breaking is not allowed in omp
    for (auto& r_node : mrModelPart.GetCommunicator().LocalMesh().Nodes()) {
        const bool is_close_enough = NodeIsCloseEnough(r_node, rThePoint, DistanceThreshold);
        if (is_close_enough) {
            local_node_found = 1;
            found_node_id = r_node.Id();
            break;
        }
    }

    CheckResults("Node", rThePoint, local_node_found);

    return found_node_id;
}

int BruteForcePointLocator::FindElement(const Point& rThePoint,
                                        Vector& rShapeFunctionValues,
                                        const double LocalCoordTol) const
{
    int found_element_id = -1; // if no element is found this will be returned
    const auto& r_elements = mrModelPart.GetCommunicator().LocalMesh().Elements();
    FindObject(r_elements, "Element",
                rThePoint, found_element_id,
                rShapeFunctionValues,
                LocalCoordTol);
    return found_element_id;
}

int BruteForcePointLocator::FindCondition(const Point& rThePoint,
                                          Vector& rShapeFunctionValues,
                                          const double LocalCoordTol) const
{
    int found_condition_id = -1; // if no condition is found this will be returned
    const auto& r_conditions = mrModelPart.GetCommunicator().LocalMesh().Conditions();
    FindObject(r_conditions, "Condition",
                rThePoint, found_condition_id,
                rShapeFunctionValues,
                LocalCoordTol);
    return found_condition_id;
}

template<typename TObjectType>
void BruteForcePointLocator::FindObject(const TObjectType& rObjects,
                                        const std::string& rObjectName,
                                        const Point& rThePoint,
                                        int& rObjectId,
                                        Vector& rShapeFunctionValues,
                                        const double LocalCoordTol) const
{
    int local_object_found = 0;
    array_1d<double, 3> local_coordinates;

    const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();
    const int global_num_objects = r_data_comm.SumAll(static_cast<int>(rObjects.size()));

    KRATOS_WARNING_IF("BruteForcePointLocator", global_num_objects == 0) << r_data_comm << "No " << rObjectName << " in ModelPart \"" << mrModelPart.Name() << std::endl;

    // note that this cannot be omp bcs breaking is not allowed in omp
    for (auto& r_object : rObjects) {
        const bool is_inside = r_object.GetGeometry().IsInside(rThePoint, local_coordinates, LocalCoordTol);
        if (is_inside) {
            local_object_found = 1;
            rObjectId = r_object.Id();
            // resizing of rShapeFunctionValues happens inside the function if required
            r_object.GetGeometry().ShapeFunctionsValues(rShapeFunctionValues, local_coordinates);
            break;
        }
    }

    CheckResults(rObjectName, rThePoint, local_object_found);
}

void BruteForcePointLocator::CheckResults(const std::string& rObjectName,
                                          const Point& rThePoint,
                                          const int LocalObjectFound) const
{
    const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();

    const int global_objects_found = r_data_comm.SumAll(LocalObjectFound);

    if (global_objects_found > 1) {
        KRATOS_WARNING_IF_ALL_RANKS("BruteForcePointLocator", LocalObjectFound == 1) << r_data_comm << "More than one " << rObjectName << " found for Point:" << rThePoint << std::endl;
    } else if (global_objects_found == 0) {
        KRATOS_WARNING("BruteForcePointLocator") << r_data_comm << "No " << rObjectName << " found for Point: " << rThePoint << std::endl;
    }
}

bool BruteForcePointLocator::NodeIsCloseEnough(const Node<3>& rNode,
                                               const Point& rThePoint,
                                               const double DistanceThreshold) const
{
    const double distance = std::sqrt( std::pow(rNode.X0() - rThePoint.X(),2)
                                     + std::pow(rNode.Y0() - rThePoint.Y(),2)
                                     + std::pow(rNode.Z0() - rThePoint.Z(),2) );

    return (distance < DistanceThreshold);
}

}  // namespace Kratos.


