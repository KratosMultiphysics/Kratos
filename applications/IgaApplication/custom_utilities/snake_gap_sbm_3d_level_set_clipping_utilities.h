//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include <array>
#include <cstddef>
#include <string>

// Project includes
#include "custom_processes/snake_gap_sbm_process.h"

namespace Kratos
{

/**
 * @brief Experimental linear level-set clipping workflow for 3D Gap-SBM.
 *
 * This utility intentionally does not create any type1/type2/type3 geometry.
 * It clips the external knot-span hexahedra against a piecewise-linear skin and
 * emits physical quadrature-point geometries for the existing Gap-SBM elements
 * and conditions.
 */
class KRATOS_API(IGA_APPLICATION) SnakeGapSbm3DLevelSetClippingUtilities
{
public:
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodePointerType = NodeType::Pointer;
    using GeometryPointerType = GeometryType::Pointer;
    using NeighbourGeometriesVectorType = std::vector<GeometryPointerType>;

    struct Settings
    {
        bool KeepNegativePhi = true;
        bool CreateInterfaceConditions = true;
        std::size_t VolumeIntegrationOrder = 2;
        std::size_t BoundaryIntegrationOrder = 2;
        std::size_t SubdivisionsPerSpan = 1;
        int EchoLevel = 0;
        std::string VolumeElementName = "GapSbmSolidElementVolumetric";
        std::string BoundaryConditionName = "GapSbmSolidCondition";
        std::string InterfaceConditionName = "GapSbmSolidInterfaceCondition";
    };

    struct Counters
    {
        std::size_t NumberOfProcessedHexCells = 0;
        std::size_t NumberOfGeneratedTetrahedraBeforeClipping = 0;
        std::size_t NumberOfClippedVolumeSubcells = 0;
        std::size_t NumberOfCutBoundaryFacets = 0;
        std::size_t NumberOfInterfaceFaces = 0;
        std::size_t NumberOfNonManifoldFaces = 0;
        std::size_t NumberOfQuadraturePointElementsCreated = 0;
        std::size_t NumberOfVolumeAuditFailures = 0;
        std::size_t NumberOfDistinctNeighbourGeometries = 0;
        std::size_t MinimumClippedTetrahedraPerNeighbour = 0;
        std::size_t MaximumClippedTetrahedraPerNeighbour = 0;
        double TotalClippedPolyhedronVolume = 0.0;
        double TotalGeneratedTetraVolume = 0.0;
        double TotalInterfaceArea = 0.0;
    };

    Counters Create(
        SnakeGapSbmProcess::Gap3DData& rGapUtilities,
        ModelPart& rRootModelPart,
        ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart,
        const SnakeGapSbmProcess::Gap3DData::KnotSpanGridInfo& rGridInfo,
        const std::set<SnakeGapSbmProcess::Gap3DData::SpanKey3D>& rActiveBackgroundSpans,
        const Vector& rKnotSpanSizes,
        ModelPart& rGapElementsModelPart,
        ModelPart& rGapConditionsModelPart,
        ModelPart& rGapInterfacesModelPart,
        Properties::Pointer pProperties,
        const Settings& rSettings) const;
};

} // namespace Kratos
