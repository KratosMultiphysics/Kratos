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

template<class TNodeType = Node>
class KRATOS_API(KRATOS_CORE) BrepSbmUtilities
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(BrepSbmUtilities);
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    using IntegrationPointType = IntegrationPoint<3>;
    using IntegrationPointsArrayType = std::vector<IntegrationPointType>;

    using GeometryType = Geometry<TNodeType>;
    using GeometryPointerType = typename GeometryType::Pointer;
    using GeometrySurrogateArrayType = DenseVector<GeometryPointerType>;


    //template<class TBrepLoopType, class TPointType>
    static void CreateBrepSurfaceSbmIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries,
        GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries,
        IntegrationInfo& rIntegrationInfo);

    
    static void CreateBrepVolumeSbmIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const std::vector<double>& rSpansW,
        GeometrySurrogateArrayType& rOuterLoops,
        GeometrySurrogateArrayType& rInnerLoops,
        IntegrationInfo& rIntegrationInfo);

private:

    static int FindKnotSpans1D(
        const std::vector<double>& rSpans, const double coord);

    ///@}
};
///@} // Kratos Classes
} // namespace Kratos.
