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
#include "utilities/reduction_utilities.h"
#include "calculate_wave_height_utility.h"


namespace Kratos
{

CalculateWaveHeightUtility::CalculateWaveHeightUtility(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart)
{
    Parameters default_parameters(R"({
        "model_part_name"  : ""
        "mean_water_level" : 0.0,
        "search_tolerance" : 1.0
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    const array_1d<double,3> gravity = mrModelPart.GetProcessInfo()[GRAVITY];
    mDirection = -gravity / norm_2(gravity);
    mMeanWaterLevel = ThisParameters["mean_water_level"].GetBool();
    mSearchTolerance = ThisParameters["search_tolerance"].GetDouble();
}

double CalculateWaveHeightUtility::Calculate(const array_1d<double,3>& rCoordinates) const
{
    KRATOS_TRY

    using MultipleReduction = CombinedReduction<SumReduction<double>, SumReduction<double>>; 

    double counter = 0.0;
    double wave_height = 0.0;

    std::tie(counter, height) = block_for_each<MultipleReduction>(
        mrModelPart.Nodes(), [&](NodeType& rNode)
    {
        double local_count = 0.0;
        double local_wave_height = 0.0;

        if (rNode->IsNot(ISOLATED) && rNode->IsNot(RIGID) && rNode->Is(FREE_SURFACE))
        {
            const auto vector_distance = MathUtils<double>::CrossProduct(mDirection - rCoordinates, rNode);
            const double distance = norm_2(vector_distance);

            if (distance < mSearchTolerance)
            {
                local_wave_height = inner_prod(mDirection, rNode) - mHeightReference;
                local_count = 1.0;
            }
        }
        double to_sum = data_vector[i];
        double to_max = data_vector[i];
        return std::make_tuple(local_count, local_wave_height);
    });

    wave_height /= counter;
    return wave_height;

    KRATOS_CATCH("");
}

}  // namespace Kratos.
