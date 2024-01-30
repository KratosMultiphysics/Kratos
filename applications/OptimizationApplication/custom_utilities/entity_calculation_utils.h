//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) EntityCalculationUtils {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = Geometry<Node>;

    using IntegrationMethodType = GeometryData::IntegrationMethod;

    ///@}
    ///@name Static operations
    ///@{

    static Geometry<Node>::Pointer CreateSolidGeometry(const GeometryType& rSurfaceGeometry);

    static void CalculateElementGaussPointData(
        Vector& rGaussWeights,
        Matrix& rShapeFunctionValues,
        const GeometryType& rGeometry,
        const IntegrationMethodType& rIntegrationMethod);

    template<class TMatrixType>
    static void CalculateSurfaceElementShapeDerivatives(
        TMatrixType& rOutput,
        GeometryType& rSolidGeometry,
        const GeometryType& rSurfaceGeometry,
        const IntegrationMethodType& rIntegrationMethod,
        const IndexType PointNumber);

    ///@}
};

///@}
} // namespace Kratos