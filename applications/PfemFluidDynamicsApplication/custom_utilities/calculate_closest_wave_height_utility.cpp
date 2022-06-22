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
#include "calculate_closest_wave_height_utility.h"

namespace Kratos
{

    CalculateClosestWaveHeightUtility::CalculateClosestWaveHeightUtility(
        ModelPart &rThisModelPart,
        Parameters ThisParameters) : mrModelPart(rThisModelPart)
    {
        Parameters default_parameters(R"({
        "model_part_name"        : "",
        "mean_water_level"       : 0.0
    })");

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        const array_1d<double, 3> gravity = mrModelPart.GetProcessInfo()[GRAVITY];
        mDirection = -gravity / norm_2(gravity);
        mMeanWaterLevel = ThisParameters["mean_water_level"].GetDouble();
        double mean_elem_size = block_for_each<SumReduction<double>>(
            mrModelPart.Elements(), [](Element &rElement)
            { return rElement.GetGeometry().Length(); });
        mean_elem_size /= mrModelPart.NumberOfElements();
        mSearchRadius = mean_elem_size ;
    }

    double CalculateClosestWaveHeightUtility::Calculate(const array_1d<double, 3> &rCoordinates) const
    {
        KRATOS_TRY

        const double safety_factor = 10; // this safety factor is to allow to find the closest node to the gauge also when mesh refinement is used 
        double search_radius = safety_factor * mSearchRadius;
        double wave_height = 0.0;
        double min_distance = search_radius;
        const auto it_node_begin = mrModelPart.NodesBegin();

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it_node = it_node_begin + i;
            if (it_node->IsNot(ISOLATED) && it_node->IsNot(RIGID) && it_node->Is(FREE_SURFACE))
            {
                const array_1d<double, 3> relative_position = it_node->Coordinates() - rCoordinates;
                const array_1d<double, 3> horizontal_position = MathUtils<double>::CrossProduct(mDirection, relative_position);
                const double distance = norm_2(horizontal_position);
                if (distance < search_radius && distance < min_distance)
                {
                    wave_height = inner_prod(mDirection, *it_node) - mMeanWaterLevel;
                    min_distance = distance;
                }
            }
        }
        return wave_height;

        KRATOS_CATCH("");
    }

} // namespace Kratos.
