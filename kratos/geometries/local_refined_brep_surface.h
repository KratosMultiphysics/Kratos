//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Wataru Fukuda
//

#pragma once

// Project includes
#include "geometries/brep_surface.h"

namespace Kratos
{

/**
 * @class LocalRefinedBrepSurface
 * @brief BrepSurface augmented with a locally-refined integration domain.
 *
 * The template parameter TLocalRefinedSurfaceType is the concrete locally-refined
 * surface class (THBSurfaceGeometry, LRSurfaceGeometry, TSurfaceGeometry, …).
 * It only needs to implement CreateIntegrationPoints and
 * CreateQuadraturePointGeometries — the same interface as any Geometry.
 *
 * Integration is fully delegated to mpLocalRefinedSurface.  The underlying
 * NurbsSurface (stored in the BrepSurface base) is used only for topology:
 * trimming-curve evaluation and geometry-part access.
 */
template<
    class TContainerPointType,
    class TLocalRefinedSurfaceType,
    bool  TShiftedBoundary = false,
    class TContainerPointEmbeddedType = TContainerPointType>
class LocalRefinedBrepSurface
    : public BrepSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LocalRefinedBrepSurface);

    using BaseType = BrepSurface<TContainerPointType, TShiftedBoundary, TContainerPointEmbeddedType>;

    using NurbsSurfaceType                = typename BaseType::NurbsSurfaceType;
    using BrepCurveOnSurfaceLoopArrayType = typename BaseType::BrepCurveOnSurfaceLoopArrayType;
    using LocalRefinedSurfaceType         = TLocalRefinedSurfaceType;

    using GeometriesArrayType             = typename BaseType::GeometriesArrayType;
    using IndexType                       = typename BaseType::IndexType;
    using IntegrationPointsArrayType      = typename BaseType::IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for untrimmed locally-refined patch.
    LocalRefinedBrepSurface(
        typename NurbsSurfaceType::Pointer          pNurbsSurface,
        typename LocalRefinedSurfaceType::Pointer   pLocalRefinedSurface)
        : BaseType(pNurbsSurface)
        , mpLocalRefinedSurface(pLocalRefinedSurface)
    {
    }

    /// Constructor for trimmed locally-refined patch.
    LocalRefinedBrepSurface(
        typename NurbsSurfaceType::Pointer          pNurbsSurface,
        typename LocalRefinedSurfaceType::Pointer   pLocalRefinedSurface,
        BrepCurveOnSurfaceLoopArrayType&            rOuterLoopArray,
        BrepCurveOnSurfaceLoopArrayType&            rInnerLoopArray)
        : BaseType(pNurbsSurface, rOuterLoopArray, rInnerLoopArray)
        , mpLocalRefinedSurface(pLocalRefinedSurface)
    {
    }

    /**
     * @brief Construct from an existing BrepSurface plus a locally-refined surface.
     *
     * Copies all topology (NurbsSurface pointer, trimming loops) from rBase via
     * BrepSurface's copy constructor, then stores pLocalRefinedSurface.
     * This is the primary factory path used in ReadLocalRefinement.
     */
    LocalRefinedBrepSurface(
        const BaseType&                             rBase,
        typename LocalRefinedSurfaceType::Pointer   pLocalRefinedSurface)
        : BaseType(rBase)
        , mpLocalRefinedSurface(pLocalRefinedSurface)
    {
    }

    /// Copy constructor.
    LocalRefinedBrepSurface(const LocalRefinedBrepSurface& rOther)
        : BaseType(rOther)
        , mpLocalRefinedSurface(rOther.mpLocalRefinedSurface)
    {
    }

    /// Destructor.
    ~LocalRefinedBrepSurface() override = default;

    ///@}
    ///@name Integration
    ///@{

    /**
     * @brief Delegates integration-point generation to the locally-refined surface.
     *
     * The locally-refined surface knows which cells are active at each level
     * and generates integration points accordingly, bypassing the NURBS span
     * decomposition of the base BrepSurface.
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        mpLocalRefinedSurface->CreateIntegrationPoints(rIntegrationPoints, rIntegrationInfo);
    }

    /**
     * @brief Builds quadrature-point geometries from the locally-refined surface.
     *
     * Shape functions are evaluated at the locally-active level for each
     * integration point.  After construction the geometry parent is set
     * so that the quadrature-point geometries point back to this BrepSurface.
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        mpLocalRefinedSurface->CreateQuadraturePointGeometries(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            rIntegrationPoints,
            rIntegrationInfo);

        for (IndexType i = 0; i < rResultGeometries.size(); ++i)
            rResultGeometries(i)->SetGeometryParent(this);
    }

    ///@}
    ///@name Access
    ///@{

    /// Returns the locally-refined surface pointer.
    typename LocalRefinedSurfaceType::Pointer pGetLocalRefinedSurface()
    {
        return mpLocalRefinedSurface;
    }

    /// Returns the locally-refined surface pointer (const).
    typename LocalRefinedSurfaceType::Pointer pGetLocalRefinedSurface() const
    {
        return mpLocalRefinedSurface;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    typename LocalRefinedSurfaceType::Pointer mpLocalRefinedSurface;

    ///@}
};

} // namespace Kratos
