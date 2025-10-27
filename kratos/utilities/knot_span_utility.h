//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:
//

#if !defined(KRATOS_KNOT_SPAN_UTILITY_H_INCLUDED)
#define KRATOS_KNOT_SPAN_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "geometries/knot_span_geometry.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    /// A Class for the creation of integration points
    template<class TPointType>
    class CreateKnotSpanUtility
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef Geometry<TPointType> GeometryType;
        typedef typename Geometry<TPointType>::Pointer GeometryPointerType;

        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef PointerVector<TPointType> PointsArrayType;

        ///@}
        ///@name Operations
        ///@{

        static GeometryPointerType CreateKnotSpanCurveOnSurface(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            double LocalTangentU,
            double LocalTangentV,
            GeometryType* pGeometryParent)
        {
            // return Kratos::make_shared<
            //     KnotSpanCurveOnSurfaceGeometry<TPointType>>(
            //         rPoints,
            //         rShapeFunctionContainer,
            //         LocalTangentU,
            //         LocalTangentV,
            //         pGeometryParent);
        }

        static GeometryPointerType CreateKnotSpanCurveOnSurface(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            double LocalTangentU,
            double LocalTangentV)
        {
            // return Kratos::make_shared<
            //     KnotSpanCurveOnSurfaceGeometry<TPointType>>(
            //         rPoints,
            //         rShapeFunctionContainer,
            //         LocalTangentU,
            //         LocalTangentV);
        }

        static GeometryPointerType CreateKnotSpanSurfaceInVolume(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            Matrix LocalTangentMatrix,
            GeometryType* pGeometryParent)
        {
            // return Kratos::make_shared<
            //     KnotSpanSurfaceInVolumeGeometry<TPointType>>(
            //         rPoints,
            //         rShapeFunctionContainer,
            //         LocalTangentMatrix,
            //         pGeometryParent);
        }

        static GeometryPointerType CreateKnotSpanSurfaceInVolume(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            Matrix LocalTangentMatrix)
        {
            // return Kratos::make_shared<
            //     KnotSpanSurfaceInVolumeGeometry<TPointType>>(
            //         rPoints,
            //         rShapeFunctionContainer,
            //         LocalTangentMatrix );
        }

        static GeometryPointerType CreateKnotSpan(
            SizeType WorkingSpaceDimension,
            SizeType LocalSpaceDimension,
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            GeometryType* pGeometryParent,
            std::size_t PointsInU,
            std::size_t PointsInV,
            double KnotSpanIntervalUBegin,
            double KnotSpanIntervalUEnd,
            double KnotSpanIntervalVBegin,
            double KnotSpanIntervalVEnd)
        {
            // if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                    KnotSpanGeometry<TPointType, 1>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent,
                        PointsInU,
                        PointsInV,
                        KnotSpanIntervalUBegin,
                        KnotSpanIntervalUEnd,
                        KnotSpanIntervalVBegin,
                        KnotSpanIntervalVEnd);
            // else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
            //     return Kratos::make_shared<
            //         KnotSpanGeometry<TPointType, 2, 1>>(
            //             rPoints,
            //             rShapeFunctionContainer,
            //             pGeometryParent);
            // else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 1)
            //     return Kratos::make_shared<
            //         KnotSpanGeometry<TPointType, 3, 1>>(
            //             rPoints,
            //             rShapeFunctionContainer,
            //             pGeometryParent);
            // else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
            //     return Kratos::make_shared<
            //         KnotSpanGeometry<TPointType, 2>>(
            //             rPoints,
            //             rShapeFunctionContainer,
            //             pGeometryParent);
            // else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
            //     return Kratos::make_shared<
            //         KnotSpanGeometry<TPointType, 3, 2>>(
            //             rPoints,
            //             rShapeFunctionContainer,
            //             pGeometryParent);
            // else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
            //     return Kratos::make_shared<
            //         KnotSpanGeometry<TPointType, 3>>(
            //             rPoints,
            //             rShapeFunctionContainer,
            //             pGeometryParent);
            // else{
            //     KRATOS_ERROR << "Working/Local space dimension combinations are "
            //         << "not provided for KnotSpanGeometry. WorkingSpaceDimension: "
            //         << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
            //         <<  std::endl;
            // }
        }

    };
    ///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_KNOT_SPAN_UTILITY_H_INCLUDED defined
