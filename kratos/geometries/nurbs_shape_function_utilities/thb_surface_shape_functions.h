//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED

// Project includes
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"

namespace Kratos {
class THBSurfaceShapeFunction
{
public:

    explicit THBSurfaceShapeFunction(
        IndexType NumberOfShapeFunctionDerivatives);

    void ComputeTHBShapeFunctionValues(
        const THBSurfaceGeometry& rGeometry,
        double Xi,
        double Eta);

    const std::vector<IndexType>&
    ControlPointIndices() const;

    SizeType NumberOfNonzeroControlPoints() const;

    double operator()(
        IndexType BasisFunctionIndex,
        IndexType DerivativeIndex) const;

private:

    void FindActiveBasisFunctions(
        const THBSurfaceGeometry& rGeometry,
        double Xi,
        double Eta);

    void EvaluateTensorProductBasisFunctions(
        const THBSurfaceGeometry& rGeometry,
        double Xi,
        double Eta);

    void ApplyTruncation(
        const THBSurfaceGeometry& rGeometry);

private:

    std::vector<IndexType> mControlPointIndices;

    Matrix mShapeFunctionData;
};
} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED defined 
