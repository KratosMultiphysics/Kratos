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

    std::pair<double,int> closest_object = {tolerance*tolerance, rObjectId};

    closest_object = block_for_each<MinDistanceReduction>(rObjects, [&](auto& r_object){

        std::vector<array_1d<double, 3>> mp_coord = { ZeroVector(3) };
        r_object.CalculateOnIntegrationPoints(MP_COORD, mp_coord, mrModelPart.GetProcessInfo());

        const double distance = std::pow(mp_coord[0][0] - rThePoint.X(),2) +
                                std::pow(mp_coord[0][1] - rThePoint.Y(),2) +
                                std::pow(mp_coord[0][2] - rThePoint.Z(),2);

        if (distance < closest_object.first) {
            closest_object = {distance, r_object.Id()};
        }

        return closest_object;
    });

    rObjectId = closest_object.second;

    KRATOS_WARNING_IF("BruteForceMaterialPointLocator", rObjectId == -1)
        << "No " << rObjectName << " close to point: " << rThePoint << std::endl;
}

}  // namespace Kratos.
