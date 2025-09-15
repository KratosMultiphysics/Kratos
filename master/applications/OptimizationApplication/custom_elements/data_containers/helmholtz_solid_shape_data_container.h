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
#include <cmath>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/constitutive_law.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_utilities/entity_calculation_utils.h"
#include "helmholtz_variable_data.h"
#include "optimization_application_variables.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class HelmholtzSolidShapeDataContainer
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = Geometry<Node>;

    static constexpr IndexType NumberOfNodes = TNumNodes;

    static constexpr IndexType NumberOfVariables = TDim;

    static constexpr auto TargetVariablesList = HelmholtzVariableData<NumberOfVariables>::TargetVariablesList;

    static constexpr auto SourceVariablesList = HelmholtzVariableData<NumberOfVariables>::SourceVariablesList;

    ///@}
    ///@name Public classes
    ///@{

    class ConstantDataContainer
    {
    public:
        ///@name Life cycle
        ///@{

        ConstantDataContainer(
            const Element& rElement,
            const GeometryData::IntegrationMethod& rIntegrationMethod,
            const ProcessInfo& rProcessInfo)
            : mrElement(rElement),
              mrIntegrationMethod(rIntegrationMethod),
              mrProcessInfo(rProcessInfo)
        {
            Vector detJ;
            mNs = mrElement.GetGeometry().ShapeFunctionsValues(rIntegrationMethod);
            mrElement.GetGeometry().ShapeFunctionsIntegrationPointsGradients(mdNdXs, detJ, mrIntegrationMethod);
        }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        const Element& mrElement;

        const GeometryData::IntegrationMethod& mrIntegrationMethod;

        const ProcessInfo& mrProcessInfo;

        GeometryType::ShapeFunctionsGradientsType mdNdXs;

        Matrix mNs;

        ///@}
        ///@name Friend classes
        ///@{

        friend class HelmholtzSolidShapeDataContainer;

        ///@}
    };

    ///@}
    ///@name Life cycle
    ///@{

    HelmholtzSolidShapeDataContainer(const Geometry<Node>& rGeometry)
    {
    }

    void AddStiffnessGaussPointContributions(
        Matrix& rStiffnessMatrix,
        const double W,
        const IndexType IntegrationPoint,
        const ConstantDataContainer& rConstantData) const
    {
        const Matrix& rdNdX = rConstantData.mdNdXs[IntegrationPoint];
        BoundedMatrix<double, 6, TNumNodes * 3> B;
        B.clear();

        IndexType local_index = 0;
        for (IndexType i = 0; i < TNumNodes; ++i) {
            B(0, local_index + 0) = rdNdX(i, 0);
            B(1, local_index + 1) = rdNdX(i, 1);
            B(2, local_index + 2) = rdNdX(i, 2);
            B(3, local_index + 0) = rdNdX(i, 1);
            B(3, local_index + 1) = rdNdX(i, 0);
            B(4, local_index + 1) = rdNdX(i, 2);
            B(4, local_index + 2) = rdNdX(i, 1);
            B(5, local_index + 0) = rdNdX(i, 2);
            B(5, local_index + 2) = rdNdX(i, 0);
            local_index += 3;
        }

        const auto& r_element = rConstantData.mrElement;
        const auto& r_geometry = r_element.GetGeometry();
        const auto& r_properties = r_element.GetProperties();

        ConstitutiveLaw::Parameters cl_params(r_geometry, r_properties, rConstantData.mrProcessInfo);
        cl_params.SetShapeFunctionsValues(row(rConstantData.mNs, IntegrationPoint));

        Matrix constitutive_matrix;
        constitutive_matrix = r_properties.GetValue(CONSTITUTIVE_LAW)->CalculateValue(cl_params, CONSTITUTIVE_MATRIX, constitutive_matrix);

        noalias(rStiffnessMatrix) += prod(trans(B), W * Matrix(prod(constitutive_matrix, B)));
    }

    ///@}
};

} // namespace Kratos