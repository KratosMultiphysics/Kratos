//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Irazabal
//                   DamApplication contributors
//

#pragma once

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/spatial_containers_configure.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

/// @brief Transfers nodal results between two model parts.
/// @details The element containing each destination node is located with the
/// Kratos Core BinBasedFastPointLocator. The Dam-specific part is limited to
/// the adaptive search tolerance and the shape-function interpolation with
/// component-wise clamping. The current configuration is used for both the
/// search and the interpolation, consistently with the Core locator.
/// @tparam TDim The spatial dimension (2 or 3).
template <SizeType TDim>
class MappingVariablesUtilities
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MappingVariablesUtilities);

    using GeometryType = Geometry<Node>;

    /// @brief Spatial configure using bounding-box overlap.
    /// @details The default Core configure tests cell overlap through
    /// Geometry::HasIntersection. Some Dam elements (for example the
    /// interface elements) use geometries that do not implement that method,
    /// so the original Dam cell search used bounding boxes. Keeping that
    /// criterion also avoids selecting cells through geometry intersections.
    class SpatialSearchConfigure : public SpatialContainersConfigure<TDim>
    {
    public:
        using BaseType = SpatialContainersConfigure<TDim>;
        using typename BaseType::PointerType;
        using typename BaseType::PointType;

        static inline bool IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint)
        {
            PointType low_point;
            PointType high_point;
            BaseType::CalculateBoundingBox(rObject, low_point, high_point);

            for (SizeType i = 0; i < TDim; ++i) {
                if (high_point[i] < rLowPoint[i] || low_point[i] > rHighPoint[i]) {
                    return false;
                }
            }
            return true;
        }

        static inline bool IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint,
            const double Radius)
        {
            return IntersectionBox(rObject, rLowPoint, rHighPoint);
        }
    };

    using LocatorType = BinBasedFastPointLocator<TDim, SpatialSearchConfigure>;

    /// @brief Maps the thermal nodal variables.
    void MappingThermalModelParts(
        ModelPart& rModelPartOld,
        ModelPart& rModelPartNew,
        const bool AddTemperature,
        const bool AddReferenceTemperature)
    {
        LocatorType locator(rModelPartOld);
        locator.UpdateSearchDatabase();

        MapNodes(
            rModelPartNew,
            locator,
            [&](Node& rNode, GeometryType& rGeometry, const Vector& rShapeFunctions) {
                if (AddTemperature) {
                    rNode.FastGetSolutionStepValue(TEMPERATURE) =
                        InterpolateAndClamp(
                            rGeometry, rShapeFunctions,
                            [](const Node& rOldNode) {
                                return rOldNode.FastGetSolutionStepValue(TEMPERATURE);
                            });
                }
                if (AddReferenceTemperature) {
                    rNode.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) =
                        InterpolateAndClamp(
                            rGeometry, rShapeFunctions,
                            [](const Node& rOldNode) {
                                return rOldNode.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE);
                            });
                }
            },
            [](Node&) {});
    }

    /// @brief Maps the mechanical nodal variables.
    void MappingMechanicalModelParts(
        ModelPart& rModelPartOld,
        ModelPart& rModelPartNew,
        const bool AddDisplacement,
        const bool AddStress)
    {
        LocatorType locator(rModelPartOld);
        locator.UpdateSearchDatabase();

        MapNodes(
            rModelPartNew,
            locator,
            [&](Node& rNode, GeometryType& rGeometry, const Vector& rShapeFunctions) {
                if (AddDisplacement) {
                    rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =
                        InterpolateAndClamp(
                            rGeometry, rShapeFunctions,
                            [](const Node& rOldNode) {
                                return rOldNode.FastGetSolutionStepValue(DISPLACEMENT_X);
                            });
                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =
                        InterpolateAndClamp(
                            rGeometry, rShapeFunctions,
                            [](const Node& rOldNode) {
                                return rOldNode.FastGetSolutionStepValue(DISPLACEMENT_Y);
                            });
                    if constexpr (TDim == 3) {
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =
                            InterpolateAndClamp(
                                rGeometry, rShapeFunctions,
                                [](const Node& rOldNode) {
                                    return rOldNode.FastGetSolutionStepValue(DISPLACEMENT_Z);
                                });
                    }
                }

                if (AddStress) {
                    Matrix nodal_stress = ZeroMatrix(TDim, TDim);
                    for (SizeType i = 0; i < TDim; ++i) {
                        for (SizeType j = 0; j < TDim; ++j) {
                            nodal_stress(i, j) = InterpolateAndClamp(
                                rGeometry, rShapeFunctions,
                                [i, j](const Node& rOldNode) {
                                    return rOldNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(i, j);
                                });
                        }
                    }
                    rNode.FastGetSolutionStepValue(INITIAL_NODAL_CAUCHY_STRESS_TENSOR) = nodal_stress;
                    rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = nodal_stress;
                }
            },
            [&](Node& rNode) {
                if (AddStress) {
                    rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = ZeroMatrix(TDim, TDim);
                    rNode.FastGetSolutionStepValue(INITIAL_NODAL_CAUCHY_STRESS_TENSOR) = ZeroMatrix(TDim, TDim);
                }
            });
    }

private:
    static constexpr SizeType mMaxNumberOfResults = 10000;

    /// @brief Locates the element containing a point.
    /// @details The search is first attempted with the geometry tolerance,
    /// then with the increased tolerances historically used by the Dam
    /// mapping utility. This keeps nodes close to the old mesh boundary
    /// mapped instead of dropping them.
    template <class TLocator>
    static bool FindElement(
        TLocator& rLocator,
        const array_1d<double, 3>& rCoordinates,
        Vector& rShapeFunctions,
        Element::Pointer& pElement,
        std::vector<Element::Pointer>& rResults)
    {
        if (rLocator.FindPointOnMesh(
                rCoordinates,
                rShapeFunctions,
                pElement,
                rResults.begin(),
                mMaxNumberOfResults,
                std::numeric_limits<double>::epsilon())) {
            return true;
        }

        for (int i = 0; i < 12; ++i) {
            for (int j = 1; j < 10; ++j) {
                const double tolerance = j * std::pow(10.0, i - 12);
                if (rLocator.FindPointOnMesh(
                        rCoordinates,
                        rShapeFunctions,
                        pElement,
                        rResults.begin(),
                        mMaxNumberOfResults,
                        tolerance)) {
                    return true;
                }
            }
        }

        pElement = nullptr;
        return false;
    }

    /// @brief Interpolates a scalar nodal value and clamps it to the local
    /// nodal range.
    template <class TValueGetter>
    static double InterpolateAndClamp(
        const GeometryType& rGeometry,
        const Vector& rShapeFunctions,
        TValueGetter GetValue)
    {
        const SizeType number_of_points = rGeometry.PointsNumber();
        Vector nodal_values(number_of_points);

        double minimum = std::numeric_limits<double>::max();
        double maximum = std::numeric_limits<double>::lowest();
        for (SizeType i = 0; i < number_of_points; ++i) {
            nodal_values[i] = GetValue(rGeometry[i]);
            minimum = std::min(minimum, nodal_values[i]);
            maximum = std::max(maximum, nodal_values[i]);
        }

        return std::clamp(inner_prod(rShapeFunctions, nodal_values), minimum, maximum);
    }

    /// @brief Locates and maps all destination nodes.
    template <class TLocator, class TFoundFunctor, class TMissingFunctor>
    static void MapNodes(
        ModelPart& rModelPartNew,
        TLocator& rLocator,
        TFoundFunctor&& rOnFound,
        TMissingFunctor&& rOnMissing)
    {
        const int number_of_nodes = static_cast<int>(rModelPartNew.Nodes().size());
        auto it_node_begin = rModelPartNew.NodesBegin();

        #pragma omp parallel
        {
            std::vector<Element::Pointer> results(mMaxNumberOfResults);
            Vector shape_functions;
            Element::Pointer p_element;

            #pragma omp for
            for (int i = 0; i < number_of_nodes; ++i) {
                auto it_node = it_node_begin + i;
                const auto& r_coordinates = it_node->Coordinates();

                if (FindElement(rLocator, r_coordinates, shape_functions, p_element, results)) {
                    rOnFound(*it_node, p_element->GetGeometry(), shape_functions);
                } else {
                    KRATOS_WARNING("MappingVariablesUtilities")
                        << "None of the old elements contains node "
                        << it_node->Id() << std::endl;
                    rOnMissing(*it_node);
                }
            }
        }
    }
};

using MappingVariables2DUtilities = MappingVariablesUtilities<2>;
using MappingVariables3DUtilities = MappingVariablesUtilities<3>;

} // namespace Kratos
