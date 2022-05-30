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
        "relative_search_radius" : 1.0,
        "search_tolerance"       : 1e-6
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    const array_1d<double,3> gravity = mrModelPart.GetProcessInfo()[GRAVITY];
    mDirection = -gravity / norm_2(gravity);
    mMeanWaterLevel = ThisParameters["mean_water_level"].GetDouble();
    double search_tolerance = ThisParameters["search_tolerance"].GetDouble();
    double relative_search_radius = ThisParameters["relative_search_radius"].GetDouble();
    double mean_elem_size = block_for_each<SumReduction<double>>(
        mrModelPart.Elements(), [](Element& rElement) {
        return rElement.GetGeometry().Length();
    });
    mean_elem_size /= mrModelPart.NumberOfElements();
    mSearchRadius = relative_search_radius * mean_elem_size + search_tolerance;
}

double CalculateWaveHeightUtility::Calculate(const array_1d<double,3>& rCoordinates) const
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

            if (distance < mSearchRadius)
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

}  // namespace Kratos.
