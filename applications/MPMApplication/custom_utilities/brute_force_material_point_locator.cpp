//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo Crescenzio
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "mpm_application_variables.h"
#include "custom_utilities/brute_force_material_point_locator.h"

namespace Kratos
{

int BruteForceMaterialPointLocator::FindElement(
    const Point& rThePoint,
    const double tolerance) const
{
    int found_element_id = -1; // if no element is found this will be returned
    const auto& r_elements = mrModelPart.Elements();
    FindObject(r_elements, "Element", rThePoint, found_element_id, tolerance);
    return found_element_id;
}

int BruteForceMaterialPointLocator::FindCondition(
    const Point& rThePoint,
    const double tolerance) const
{
    int found_condition_id = -1; // if no condition is found this will be returned
    const auto& r_conditions = mrModelPart.Conditions();
    FindObject(r_conditions, "Condition", rThePoint, found_condition_id, tolerance);
    return found_condition_id;
}

template<typename TObjectType>
void BruteForceMaterialPointLocator::FindObject(
    const TObjectType& rObjects,
    const std::string& rObjectName,
    const Point& rThePoint,
    int& rObjectId,
    const double tolerance) const
{
    const int num_objects = rObjects.size();

    if (num_objects == 0) {
        KRATOS_WARNING("BruteForceMaterialPointLocator") << "No " << rObjectName
            << " in ModelPart \"" << mrModelPart.Name() << "\"" << std::endl;
        return;
    }

    auto& process_info = mrModelPart.GetProcessInfo();

    const std::tuple<double,int> closest_object{ tolerance*tolerance, -1 };
    std::vector<array_1d<double, 3>> mp_coord = { ZeroVector(3) };

    rObjectId = block_for_each<MPMMinDistanceReduction>(rObjects, mp_coord,
        [&process_info, &closest_object, &rThePoint](auto& r_object, auto& r_mp_coord){

            r_object.CalculateOnIntegrationPoints(MP_COORD, r_mp_coord, process_info);

            const double distance = std::pow(r_mp_coord[0][0] - rThePoint.X(),2) +
                                    std::pow(r_mp_coord[0][1] - rThePoint.Y(),2) +
                                    std::pow(r_mp_coord[0][2] - rThePoint.Z(),2);

            if (distance < get<0>(closest_object)) {
                return std::tuple<double,int>(distance, r_object.Id());
            }

            return closest_object;
        });

    KRATOS_WARNING_IF("BruteForceMaterialPointLocator", rObjectId == -1)
        << "No " << rObjectName << " close to point: " << rThePoint
        << " (tolerance: " << tolerance << ")." << std::endl;
}

}  // namespace Kratos.
