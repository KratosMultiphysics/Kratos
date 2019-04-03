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
    int local_nodes_found = 0;

    // note that this cannot be omp bcs breaking is not allowed in omp
    for (auto& r_node : mrModelPart.GetCommunicator().LocalMesh().Nodes())
    {
        const bool is_close_enough = NodeIsCloseEnough(r_node, rThePoint, DistanceThreshold);
        if (is_close_enough)
        {
            local_nodes_found = 1;
            found_node_id = r_node.Id();
            break;
        }
    }

    CheckResults("Node", rThePoint, local_nodes_found);

    return found_node_id;
}

int BruteForcePointLocator::FindElement(const Point& rThePoint,
                              Vector& rShapeFunctionValues) const
{
    int found_element_id = -1; // if no element is found this will be returned
    const auto& r_elements = mrModelPart.GetCommunicator().LocalMesh().Elements();
    FindObject(r_elements, "Element",
                rThePoint, found_element_id,
                rShapeFunctionValues);
    return found_element_id;
}

int BruteForcePointLocator::FindCondition(const Point& rThePoint,
                                Vector& rShapeFunctionValues) const
{
    int found_condition_id = -1; // if no condition is found this will be returned
    const auto& r_conditions = mrModelPart.GetCommunicator().LocalMesh().Conditions();
    FindObject(r_conditions, "Condition",
                rThePoint, found_condition_id,
                rShapeFunctionValues);
    return found_condition_id;
}

template<typename TObjectType>
void BruteForcePointLocator::FindObject(const TObjectType& rObjects,
                              const std::string& rObjectName,
                              const Point& rThePoint,
                              int& rObjectId,
                              Vector& rShapeFunctionValues) const
{

    const int domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    const auto& r_geom = rObjects.begin()->GetGeometry();

    KRATOS_ERROR_IF_NOT(static_cast<std::size_t>(domain_size) == r_geom.WorkingSpaceDimension())
        << "Domain size (" << domain_size << ") and WorkingSpaceDimension ("
        << r_geom.WorkingSpaceDimension() << ") of the " << rObjectName
        << " are not equal!" << std::endl;

    int local_objects_found = 0;
    array_1d<double, 3> local_coordinates;

    // note that this cannot be omp bcs breaking is not allowed in omp
    for (auto& r_object : rObjects)
    {
        const bool is_inside = r_object.GetGeometry().IsInside(rThePoint, local_coordinates);
        if (is_inside)
        {
            local_objects_found = 1;
            rObjectId = r_object.Id();
            // resizing of rShapeFunctionValues happens inside the function if required
            r_object.GetGeometry().ShapeFunctionsValues(rShapeFunctionValues, local_coordinates);
            break;
        }
    }

    CheckResults(rObjectName, rThePoint, local_objects_found);
}

void BruteForcePointLocator::CheckResults(const std::string& rObjectName,
                                const Point& rThePoint,
                                int LocalObjectsFound) const
{
    mrModelPart.GetCommunicator().SumAll(LocalObjectsFound);

    if (LocalObjectsFound > 1)
    {
        KRATOS_WARNING_IF("Point Locator", mrModelPart.GetCommunicator().MyPID() == 0)
            << "More than one " << rObjectName << " found for Point: " << rThePoint << std::endl;
        mrModelPart.GetCommunicator().Barrier();
        KRATOS_WARNING("Point Locator")
            << "    In Rank: " << mrModelPart.GetCommunicator().MyPID() << std::endl;
        mrModelPart.GetCommunicator().Barrier();
    }
    else if (LocalObjectsFound == 0)
    {
        KRATOS_WARNING_IF("Point Locator", mrModelPart.GetCommunicator().MyPID() == 0)
            << "No " << rObjectName << " found for Point: " << rThePoint << std::endl;
    }
}

bool BruteForcePointLocator::NodeIsCloseEnough(const Node<3>& rNode,
                                     const Point& rThePoint,
                                     double DistanceThreshold) const
{
    const double distance = std::sqrt( std::pow(rNode.X0() - rThePoint.X(),2)
                                     + std::pow(rNode.Y0() - rThePoint.Y(),2)
                                     + std::pow(rNode.Z0() - rThePoint.Z(),2) );

    return (distance < DistanceThreshold);
}

}  // namespace Kratos.


