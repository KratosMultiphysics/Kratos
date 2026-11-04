//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication contributors
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/binbased_fast_point_locator.h"

// Application includes
#include "custom_utilities/mapping_variables_utilities.hpp"

namespace Kratos
{

/// @brief Thin wrapper around the Core point locator for Dam initial-state mapping.
/// @details The Core locator is configured with the Dam bounding-box search
/// configuration. This is required because several supported Dam geometries do
/// not implement Geometry::HasIntersection, which the default Core configure
/// relies on. The wrapper only adds local-coordinate retrieval and a strict
/// inside check to the existing locator; it does not implement a new search.
template <SizeType TDim>
class DamInitialStateSpatialSearchUtility
{
public:
    using ConfigureType = typename MappingVariablesUtilities<TDim>::SpatialSearchConfigure;
    using LocatorType = BinBasedFastPointLocator<TDim, ConfigureType>;

    struct Result
    {
        bool Found = false;
        Vector ShapeFunctions;
        Element::Pointer pElement = nullptr;
        Vector LocalCoordinates;
    };

    explicit DamInitialStateSpatialSearchUtility(ModelPart& rSourceModelPart)
        : mLocator(rSourceModelPart)
    {
    }

    void UpdateSearchDatabase()
    {
        mLocator.UpdateSearchDatabase();
    }

    /// @brief Locates a physical point and returns shape functions and local coordinates.
    /// @details The point is accepted only if the containing geometry reports it
    /// inside with the requested local-coordinate tolerance.
    Result FindPoint(
        const array_1d<double, 3>& rCoordinates,
        const double Tolerance = 1.0e-9)
    {
        Result result;
        const bool is_found = mLocator.FindPointOnMeshSimplified(
            rCoordinates, result.ShapeFunctions, result.pElement,
            mMaxNumberOfResults, Tolerance);
        if (!is_found || result.pElement == nullptr) {
            return result;
        }

        Geometry<Node>::CoordinatesArrayType local_coordinates;
        if (!result.pElement->GetGeometry().IsInside(
                rCoordinates, local_coordinates, Tolerance)) {
            return result;
        }

        result.LocalCoordinates.resize(3, false);
        for (std::size_t i = 0; i < 3; ++i) {
            result.LocalCoordinates[i] = local_coordinates[i];
        }
        result.Found = true;
        return result;
    }

private:
    static constexpr std::size_t mMaxNumberOfResults = 10000;

    LocatorType mLocator;
};

} // namespace Kratos
