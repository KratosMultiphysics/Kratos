//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#pragma once

// Std includes
#include <list>

// System includes
#include "includes/define.h"

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"

#include "includes/node.h"

#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

template<bool TShiftedBoundary, class TNodeType = Node>
class KRATOS_API(KRATOS_CORE) BrepSbmUtilities
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(BrepSbmUtilities);
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

    using BrepCurveOnSurfacePointerType = typename BrepCurveOnSurface<PointerVector<Node>, TShiftedBoundary, PointerVector<Point>>::Pointer;

    using GeometryType = Geometry<TNodeType>;
    using GeometryPointerType = typename GeometryType::Pointer;
    using GeometrySurrogateArrayType = DenseVector<GeometryPointerType>;

    /**
     * @brief Create a Brep Surface Sbm Integration Points object
     * 
     * @param rIntegrationPoints 
     * @param rSpansU 
     * @param rSpansV 
     * @param rSurrogateOuterLoopGeometries 
     * @param rSurrogateInnerLoopGeometries 
     * @param rIntegrationInfo 
     */
    static void CreateBrepSurfaceSbmIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries,
        GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries,
        IntegrationInfo& rIntegrationInfo);

    /**
     * @brief Create a Brep Volume Sbm Integration Points object
     * 
     * @param rIntegrationPoints 
     * @param rSpansU 
     * @param rSpansV 
     * @param rSpansW 
     * @param rOuterLoops 
     * @param rInnerLoops 
     * @param rIntegrationInfo 
     */
    static void CreateBrepVolumeSbmIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const std::vector<double>& rSpansW,
        GeometrySurrogateArrayType& rOuterLoops,
        GeometrySurrogateArrayType& rInnerLoops,
        IntegrationInfo& rIntegrationInfo);

private:
    
    /**
     * @brief 
     * 
     * @param rSpans 
     * @param coord 
     * @return int 
     */
    static int FindKnotSpans1D(
        const std::vector<double>& rSpans, const double coord);

    ///@}
};
///@} // Kratos Classes
} // namespace Kratos.
