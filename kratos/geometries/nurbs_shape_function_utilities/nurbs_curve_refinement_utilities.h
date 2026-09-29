//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/nurbs_curve_geometry.h"
#include "nurbs_utilities.h"
#include "includes/node.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) NurbsCurveRefinementUtilities{
    public:
        using IndexType = std::size_t;
        using SizeType = std::size_t;
        using NodeType = Node;

        using NurbsCurveGeometryType = NurbsCurveGeometry<3, PointerVector<NodeType>>;
        using NurbsCurveGeometryPointerType = typename NurbsCurveGeometryType::Pointer;

        static void KnotRefinement(
            NurbsCurveGeometryType& rGeometry,
            std::vector<double>& rKnotsToInsert,
            PointerVector<NodeType>& rPointsRefined,
            Vector& rKnotsRefined,
            Vector& rWeightsRefined);

        static void DegreeElevation(
            NurbsCurveGeometryType& rGeometry,
            const SizeType rDegreeToElevate,
            PointerVector<NodeType>& rPointsRefined,
            Vector& rKnotsRefined,
            Vector& rWeightsRefined);

        static void SortAndFilter(
            std::vector<double>& rKnotsToInsert,
            const Vector& rKnotsOld);

        static Vector CreateFullKnotVector(
            const Vector& rKnots,
            const SizeType Degree);

        static Vector CreateReducedKnotVector(
            const Vector& rFullKnots,
            const SizeType Degree);
};

} // namespace Kratos