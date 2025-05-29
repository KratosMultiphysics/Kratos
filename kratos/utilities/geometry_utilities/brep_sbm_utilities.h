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

    /**
     * @brief Generates integration points for a BREP surface using SBM (Structured Background Mesh) logic.
     * 
     * @param rSpansU Knot spans in the U parametric direction.
     * @param rSpansV Knot spans in the V parametric direction.
     * @param rSurrogateOuterLoopGeometries Geometries approximating the outer boundary loop of the surface.
     * @param rSurrogateInnerLoopGeometries Geometries approximating the inner holes.
     * @param rIntegrationPoints Output vector of computed integration points on the surface.
     * @param rIntegrationInfo Object containing integration settings and metadata.
     */
    static void CreateBrepSurfaceSbmIntegrationPoints(
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries,
        const GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries,
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo);

    /**
     * @brief Generates integration points for a BREP volume using SBM (Structured Background Mesh) logic.
     * 
     * @param rIntegrationPoints Output vector to store the computed integration points within the volume.
     * @param rSpansU Knot spans in the U parametric direction (first parametric axis).
     * @param rSpansV Knot spans in the V parametric direction (second parametric axis).
     * @param rSpansW Knot spans in the W parametric direction (third parametric axis).
     * @param rOuterLoops Geometries representing the outer boundary surfaces of the volume.
     * @param rInnerLoops Geometries representing internal holes or voids in the volume.
     * @param rIntegrationInfo Object containing integration settings.
     */
    static void CreateBrepVolumeSbmIntegrationPoints(
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const std::vector<double>& rSpansW,
        GeometrySurrogateArrayType& rOuterLoops,
        GeometrySurrogateArrayType& rInnerLoops,
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo);

private:
    
    /**
     * @brief Finds the knot span index in a 1D parametric space that contains the given coordinate.
     * 
     * @param rSpans A sorted vector of knot span boundary values.
     * @param coord The parametric coordinate for which the knot span index is to be found.
     * @return int The index of the knot span in which the coordinate lies.
     */
    static int FindKnotSpans1D(
        const std::vector<double>& rSpans, 
        const double coord);

    ///@}
};
///@} // Kratos Classes
} // namespace Kratos.
