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
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"

// Application utilities
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

    // The average element size is an initial guess for computing the average height or to finding the closest node
    mMeanElementSize = block_for_each<SumReduction<double>>(
        mrModelPart.Elements(), [](Element& rElement) {
        return rElement.GetGeometry().Length();
    });
    mMeanElementSize /= mrModelPart.NumberOfElements();

    // The top and low bounds will be used to find the local element size
    using MinMaxReduction = CombinedReduction<MinReduction<double>,MaxReduction<double>>; 

    std::tie(mLowBound, mTopBound) = block_for_each<MinMaxReduction>(mrModelPart.Nodes(), [&](NodeType& rNode){
        double vertical_position = inner_prod(mDirection, rNode.Coordinates());
        return std::make_tuple(vertical_position, vertical_position);
    });
    double distance = mTopBound - mLowBound; // increasing the limits for safety
    mTopBound += distance;
    mLowBound -= distance;
}


double CalculateWaveHeightUtility::Calculate(const array_1d<double,3>& rCoordinates) const
{
    // First step: calculating the search radius
    double search_radius;
    if (mUseLocalElementSize) {
        double element_size = FindLocalElementSize(rCoordinates);
        search_radius = mRelativeRadius * element_size + mAbsoluteRadius;
    } else {
        search_radius = mRelativeRadius * mMeanElementSize + mAbsoluteRadius;
    }

    // Second step: calculating the wave height
    if (mUseNearestNode) {
        return CalculateNearest(rCoordinates, search_radius);
    } else {
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


double CalculateWaveHeightUtility::CalculateNearest(const array_1d<double,3>& rCoordinates, double SearchRadius) const
{
    KRATOS_TRY

    struct SurfacePoint{
        double distance;
        double mean_water_level;
        array_1d<double,3> coordinates;
        array_1d<double,3> direction;
    };

    struct CustomReducer{
        using return_type = double;
        double distance = 1e16;
        double wave_height = 0.0;

        return_type GetValue()
        {
            return wave_height;
        }

        void LocalReduce(SurfacePoint& r_point)
        {
            if (r_point.distance < this->distance) {
                this->distance = r_point.distance;
                this->wave_height = inner_prod(r_point.direction, r_point.coordinates) - r_point.mean_water_level;
            }
        }

        void ThreadSafeReduce(CustomReducer& rOther)
        {
            #pragma omp critical
            {
                if (rOther.distance < this->distance) {
                    this->distance = rOther.distance;
                    this->wave_height = rOther.wave_height;
                }
            }
        }
    };

    return block_for_each<CustomReducer>(mrModelPart.Nodes(), [&](NodeType& rNode) {
        SurfacePoint surface_point;
        if (rNode.IsNot(ISOLATED) && rNode.IsNot(RIGID) && rNode.Is(FREE_SURFACE))
        {
            const array_1d<double,3> relative_position = rNode.Coordinates() - rCoordinates;
            const array_1d<double,3> horizontal_position = MathUtils<double>::CrossProduct(mDirection, relative_position);
            surface_point.distance = norm_2(horizontal_position);
            surface_point.coordinates = rNode.Coordinates();
            surface_point.direction = mDirection;
            surface_point.mean_water_level = mMeanWaterLevel;
        }
        return surface_point;
    });

    KRATOS_CATCH("")
}


double CalculateWaveHeightUtility::FindLocalElementSize(const array_1d<double,3>& rCoordinates) const
{
    // TODO: 1. move tetrahedra-line intersection to intersection utilities
    // TODO: 2. use intersection utilities and avoid the creation of nodes, use points instead
    NodeType::Pointer low_point = Kratos::make_intrusive<NodeType>(0, rCoordinates + mDirection * (mLowBound - inner_prod(rCoordinates, mDirection)));
    NodeType::Pointer top_point = Kratos::make_intrusive<NodeType>(0, rCoordinates + mDirection * (mTopBound - inner_prod(rCoordinates, mDirection)));
    GeometryType::Pointer p_vertical_line;
    if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
        p_vertical_line = Kratos::make_shared<Line3D2<NodeType>>(low_point, top_point);
    }
    else if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2) {
        p_vertical_line = Kratos::make_shared<Line2D2<NodeType>>(low_point, top_point);
    }
    else {
        KRATOS_ERROR << "CalculateWaveHeightUtility. DOMAIN_SIZE is not defined in the model part." << std::endl;
    }
    return block_for_each<MinReduction<double>>(mrModelPart.Elements(), [&](Element& rElement){
        double local_elem_size = 1e16;
        if (rElement.GetGeometry().HasIntersection(*p_vertical_line)) {
            local_elem_size = rElement.GetGeometry().Length();
        }
        return local_elem_size;
    });
}

}  // namespace Kratos.
