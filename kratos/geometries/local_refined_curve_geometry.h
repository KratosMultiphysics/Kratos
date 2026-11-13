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
#include "geometries/geometry.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class LocalRefinedCurveGeometry : public Geometry<typename TContainerPointType::value_type>
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

    KRATOS_CLASS_POINTER_DEFINITION(LocalRefinedCurveGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    explicit LocalRefinedCurveGeometry(const PointsArrayType& rThisPoints)
        : BaseType(rThisPoints, &msGeometryData)
    {}

    ~LocalRefinedCurveGeometry() override = default;

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
    ///@name Interface for locally refined curve geometries
    ///@{

    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override = 0;

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

    BaseType* mpGeometryParent = nullptr;

    ///@}

}; // class LocalRefinedCurveGeometry

template<int TWorkingSpaceDimension, class TContainerPointType>
const GeometryData LocalRefinedCurveGeometry<TWorkingSpaceDimension, TContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TContainerPointType>
const GeometryDimension LocalRefinedCurveGeometry<TWorkingSpaceDimension, TContainerPointType>::msGeometryDimension(
    TWorkingSpaceDimension, 1);

} // namespace Kratos
