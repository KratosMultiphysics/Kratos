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
#include "point_locator.h"


namespace Kratos
{
    bool PointLocator::FindNode(const Point& rThePoint,
                                int& rNodeId,
                                double DistanceThreshold)
    {
        rNodeId = -1;
        bool is_close_enough = false;

        int global_nodes_found = 0;

        // note that this cannot be omp bcs breaking is not allowed in omp
        for (auto& r_node : mrModelPart.GetCommunicator().LocalMesh().Nodes())
        {
            is_close_enough = NodeIsCloseEnough(r_node, rThePoint, DistanceThreshold);
            if (is_close_enough)
            {
                global_nodes_found = 1;
                rNodeId = r_node.Id();
                break;
            }
        }

        CheckResults("Node", rThePoint, global_nodes_found);

        return is_close_enough;
    }

    bool PointLocator::FindElement(const Point& rThePoint,
                                   int& rObjectId,
                                   Vector& rShapeFunctionValues)
    {
        const auto& r_elements = mrModelPart.GetCommunicator().LocalMesh().Elements();
        const bool is_inside = FindObject(r_elements, "Element",
                                          rThePoint, rObjectId,
                                          rLocalCoordinates);
        return is_inside;
    }

    bool PointLocator::FindCondition(const Point& rThePoint,
                                     int& rObjectId,
                                     Vector& rShapeFunctionValues)
    {
        const auto& r_conditions = mrModelPart.GetCommunicator().LocalMesh().Conditions();
        const bool is_inside = FindObject(r_conditions, "Condition",
                                          rThePoint, rObjectId,
                                          rLocalCoordinates);
        return is_inside;
    }

    template<typename TObjectType>
    bool PointLocator::FindObject(const TObjectType& rObjects,
                                  const std::string& rObjectType,
                                  const Point& rThePoint,
                                  int& rObjectId,
                                  Vector& rShapeFunctionValues)
    {

        const int domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

        const auto& r_geom = rObjects.begin()->GetGeometry();

        KRATOS_ERROR_IF_NOT(static_cast<std::size_t>(domain_size) == r_geom.WorkingSpaceDimension())
            << "Domain size (" << domain_size << ") and WorkingSpaceDimension ("
            << r_geom.WorkingSpaceDimension() << ") of the "
            << "Elements are not equal!" << std::endl;

        rObjectId = -1;
        bool is_inside;

        int global_objects_found = 0;
        array_1d<double, 3> local_coordinates;

        // note that this cannot be omp bcs breaking is not allowed in omp
        for (auto& r_object : rObjects)
        {
            is_inside = r_object.GetGeometry().IsInside(rThePoint, local_coordinates);
            if (is_inside)
            {
                global_objects_found = 1;
                rObjectId = r_object.Id();
                // resizing of rShapeFunctionValues happens inside the function if required
                r_object.GetGeometry().ShapeFunctionsValues(rShapeFunctionValues, local_coordinates);
                break;
            }
        }

        CheckResults(rObjectType, rThePoint, global_objects_found);

        return is_inside;
    }

    void PointLocator::CheckResults(const std::string& rObjectType,
                                    const Point& rThePoint,
                                    int GlobalObjectsFound)
    {
        mrModelPart.GetCommunicator().SumAll(GlobalObjectsFound);

        if (GlobalObjectsFound > 1)
        {
            KRATOS_WARNING_IF("Point Locator", mrModelPart.GetCommunicator().MyPID() == 0)
                << "More than one " << rObjectType << " found for Point: " << rThePoint << std::endl;
            mrModelPart.GetCommunicator().Barrier();
            KRATOS_WARNING("Point Locator")
                << "    In Rank: " << mrModelPart.GetCommunicator().MyPID() << std::endl;
            mrModelPart.GetCommunicator().Barrier();
        }
        else if (GlobalObjectsFound == 0)
        {
            KRATOS_WARNING_IF("Point Locator", mrModelPart.GetCommunicator().MyPID() == 0)
                << "No " << rObjectType << " found for Point: " << rThePoint << std::endl;
        }
    }


    bool PointLocator::NodeIsCloseEnough(const Node<3>& rNode,
                                         const Point& rThePoint,
                                         double DistanceThreshold)
    {
        const double distance = std::sqrt( std::pow(rNode.X() - rThePoint.X(),2)
                                         + std::pow(rNode.Y() - rThePoint.Y(),2)
                                         + std::pow(rNode.Z() - rThePoint.Z(),2) );

        return (distance < DistanceThreshold);
    }

}  // namespace Kratos.


