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
        "model_part_name"        : "",
        "mean_water_level"       : 0.0,
        "relative_search_radius" : 2.0,
        "search_tolerance"       : 1e-6,
        "use_local_element_size" : false,
        "use_nearest_node"       : false
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    const array_1d<double,3> gravity = mrModelPart.GetProcessInfo()[GRAVITY];
    mDirection = -gravity / norm_2(gravity);
    mUseLocalElementSize = ThisParameters["use_local_element_size"].GetBool();
    mUseNearestNode = ThisParameters["use_nearest_node"].GetBool();
    mMeanWaterLevel = ThisParameters["mean_water_level"].GetDouble();
    mAbsoluteRadius = ThisParameters["search_tolerance"].GetDouble();
    mRelativeRadius = ThisParameters["relative_search_radius"].GetDouble();
    double elem_size_sum = block_for_each<SumReduction<double>>(
        mrModelPart.Elements(), [](Element& rElement) {
        return rElement.GetGeometry().Length();
    });
    mMeanElementSize = elem_size_sum / mrModelPart.NumberOfElements();
}


double CalculateWaveHeightUtility::Calculate(const array_1d<double,3>& rCoordinates) const
{
    if (mUseNearestNode) {
        return CalculateNearest(rCoordinates);
    } else {
        double search_radius;
        if (mUseLocalElementSize) {
            double element_size = FindLocalElementSize(rCoordinates);
            search_radius = mRelativeRadius * element_size + mAbsoluteRadius;
        } else {
            search_radius = mRelativeRadius * mMeanElementSize + mAbsoluteRadius;
        }
        return CalculateAverage(rCoordinates, search_radius);
    }
}


double CalculateWaveHeightUtility::CalculateAverage(const array_1d<double,3>& rCoordinates, double SearchRadius) const
{
    KRATOS_TRY

    using MultipleReduction = CombinedReduction<SumReduction<double>, SumReduction<double>>; 

    double counter = 0.0;
    double wave_height = 0.0;

    std::tie(counter, wave_height) = block_for_each<MultipleReduction>(
        mrModelPart.Nodes(), [&](NodeType& rNode)
    {
        double local_count = 0.0;
        double local_wave_height = 0.0;

        if (rNode.IsNot(ISOLATED) && rNode.IsNot(RIGID) && rNode.Is(FREE_SURFACE))
        {
            const array_1d<double,3> relative_position = rNode.Coordinates() - rCoordinates;
            const array_1d<double,3> horizontal_position = MathUtils<double>::CrossProduct(mDirection, relative_position);
            const double distance = norm_2(horizontal_position);

            if (distance < SearchRadius)
            {
                local_count = 1.0;
                local_wave_height = inner_prod(mDirection, rNode) - mMeanWaterLevel;
            }
        }
        return std::make_tuple(local_count, local_wave_height);
    });

    wave_height /= counter;
    return wave_height;

    KRATOS_CATCH("");
}


double CalculateWaveHeightUtility::CalculateNearest(const array_1d<double,3>& rCoordinates) const
{
    KRATOS_TRY

    return 0.0;

    KRATOS_CATCH("")
}


double CalculateWaveHeightUtility::FindLocalElementSize(const array_1d<double,3>& rCoordinates) const
{
    return 0.0;
}


void CalculateWaveHeightUtility::UpdateSearchDatabase()
{

}

}  // namespace Kratos.
