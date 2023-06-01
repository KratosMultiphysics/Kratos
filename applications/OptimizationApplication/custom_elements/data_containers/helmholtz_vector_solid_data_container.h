//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <array>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/entity_calculation_utils.h"
#include "optimization_application_variables.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class HelmholtzVectorSolidDataContainer {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = Geometry<Node>;

    static constexpr IndexType NumberOfNodes = TNumNodes;

    static constexpr auto TargetVariablesList = (TDim == 2) ?
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_X, &HELMHOLTZ_VECTOR_Y} :
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_X, &HELMHOLTZ_VECTOR_Y, &HELMHOLTZ_VECTOR_Z};

    static constexpr auto SourceVariablesList = (TDim == 2) ?
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_SOURCE_X, &HELMHOLTZ_VECTOR_SOURCE_Y} :
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_SOURCE_X, &HELMHOLTZ_VECTOR_SOURCE_Y, &HELMHOLTZ_VECTOR_SOURCE_Z};


    ///@}
    ///@name Public classes
    ///@{

    struct ConstantDataContainer
    {
        const GeometryType& mrGeometry;
        const GeometryData::IntegrationMethod& mrIntegrationMethod;
        GeometryType::ShapeFunctionsGradientsType mdNdXs;

        ConstantDataContainer(
            const GeometryType& rGeometry,
            const GeometryData::IntegrationMethod& rIntegrationMethod)
            : mrGeometry(rGeometry),
              mrIntegrationMethod(rIntegrationMethod)
        {
        }
    };

    ///@}
    ///@name Life cycle
    ///@{

    HelmholtzVectorSolidDataContainer(const Geometry<Node>& rGeometry)
    {
    }

    void CalculateConstants(ConstantDataContainer& rConstantData) const
    {
        Vector detJ;
        rConstantData.mrGeometry.ShapeFunctionsIntegrationPointsGradients(rConstantData.mdNdXs, detJ, rConstantData.mrIntegrationMethod);
    }

    void CalculateShapeFunctionDerivatives(
        Matrix& rdNdX,
        const IndexType IntegrationPoint,
        const ConstantDataContainer& rConstantData) const
    {
        rdNdX = rConstantData.mdNdXs[IntegrationPoint];
    }

    ///@}
};

} // namespace Kratos