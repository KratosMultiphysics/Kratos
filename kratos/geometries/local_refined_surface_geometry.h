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

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class LocalRefinedSurfaceGeometry : public Geometry<typename TContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TContainerPointType::value_type NodeType;

    typedef Geometry<NodeType> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    using BaseType::CreateQuadraturePointGeometries;
    using BaseType::pGetPoint;
    using BaseType::GetPoint;

    KRATOS_CLASS_POINTER_DEFINITION(LocalRefinedSurfaceGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of points
    explicit LocalRefinedSurfaceGeometry(const PointsArrayType& rThisPoints)
        : BaseType(rThisPoints, &msGeometryData)
    {}

    ~LocalRefinedSurfaceGeometry() override = default;

    ///@}
    ///@name Parent
    ///@{

    BaseType& GetGeometryParent(IndexType Index) const override
    {
        return *mpGeometryParent;
    }

    void SetGeometryParent(BaseType* pGeometryParent) override
    {
        mpGeometryParent = pGeometryParent;
    }

    ///@}
    ///@name Interface for locally refined geometries
    ///@{

    virtual void GetActiveBoundaryPoints(
        double LocalCoordinateU,
        double LocalCoordinateV,
        IndexType VariationU,
        IndexType VariationV,
        std::vector<typename NodeType::Pointer>& rPoints) const = 0;

    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override = 0;

    SizeType PointsNumberInDirection(IndexType LocalDirectionIndex) const override = 0;

    void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType DirectionIndex) const override = 0;

    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override = 0;

    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override = 0;

    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates) const override = 0;

    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override = 0;

    ///@}

private:
    ///@name Private Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Private Member Variables
    ///@{

    /// Pointer to the wrapping BrepSurface (or LocalRefinedBrepSurface).
    BaseType* mpGeometryParent = nullptr;

    ///@}

}; // class LocalRefinedSurfaceGeometry

template<int TWorkingSpaceDimension, class TContainerPointType>
const GeometryData LocalRefinedSurfaceGeometry<TWorkingSpaceDimension, TContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TContainerPointType>
const GeometryDimension LocalRefinedSurfaceGeometry<TWorkingSpaceDimension, TContainerPointType>::msGeometryDimension(
    TWorkingSpaceDimension, 2);

} // namespace Kratos

