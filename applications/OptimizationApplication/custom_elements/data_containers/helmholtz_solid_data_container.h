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

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"

// Application includes
#include "custom_utilities/entity_calculation_utils.h"
#include "helmholtz_variable_data.h"
#include "optimization_application_variables.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TDataDimension>
class HelmholtzSolidDataContainer
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = Geometry<Node>;

    static constexpr IndexType NumberOfNodes = TNumNodes;

    static constexpr IndexType NumberOfVariables = (TDataDimension == 1) ? 1 : TDim;

    static constexpr auto TargetVariablesList = HelmholtzVariableData<NumberOfVariables>::TargetVariablesList;

    static constexpr auto SourceVariablesList = HelmholtzVariableData<NumberOfVariables>::SourceVariablesList;

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

    HelmholtzSolidDataContainer(const Geometry<Node>& rGeometry)
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