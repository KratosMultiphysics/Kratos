//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
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

    ///@}
    ///@name Static operations
    ///@{

    static Geometry<Node>::Pointer CreateSolidGeometry(const Geometry<Node>& rSurfaceGeometry);

    template<class TMatrixType>
    static void CalculateSurfaceElementShapeDerivatives(
        TMatrixType& rOutput,
        Geometry<Node>& rSolidGeometry,
        const Geometry<Node>& rSurfaceGeometry,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        const IndexType PointNumber);

    ///@}
};

///@}
} // namespace Kratos