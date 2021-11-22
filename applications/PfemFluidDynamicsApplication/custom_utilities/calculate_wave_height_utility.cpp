//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "calculate_wave_height_utility.h"


namespace Kratos
{

CalculateWaveHeightUtility::CalculateWaveHeightUtility(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart)
{
    Parameters default_parameters(R"({
        "coordinates"      : [0.0, 0.0, 0.0],
        "mean_water_level" : 0.0,
        "search_tolerance" : 1.0
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    const array_1d<double,3> gravity = mrModelPart.GetProcessInfo()[GRAVITY];
    mDirection = gravity / norm_2(gravity);
    mCoordinates = ThisParameters["coordinates"].GetBool();
    mMeanWaterLevel = ThisParameters["mean_water_level"].GetBool();
    mSearchTolerance = ThisParameters["search_tolerance"].GetDouble();
}

double CalculateWaveHeightUtility::Execute() const
{
    KRATOS_TRY

    double counter = 0.0;
    double heightSum = 0.0;

    for (NodeType& rNode : mrModelPart.Nodes())
    {
        if (rNode->IsNot(ISOLATED) && rNode->IsNot(RIGID) && rNode->Is(FREE_SURFACE))
        {
            const auto vector_distance = MathUtils<double>::CrossProduct(mDirection, rNode);
            const double distance = norm_2(vector_distance);
            if (distance < mSearchTolerance)
            {
                const double wave_height = inner_prod(mDirection, rNode) - mHeightReference;
                counter += 1.0;
                heightSum += wave_height;
            }
        }
    }

    const double height = heightSum / counter;
    return height;

    KRATOS_CATCH("");
}

}  // namespace Kratos.
