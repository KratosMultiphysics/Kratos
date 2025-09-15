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
#include "includes/element.h"
#include "includes/process_info.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/entity_calculation_utils.h"
#include "helmholtz_variable_data.h"
#include "optimization_application_variables.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TDataDimension>
class HelmholtzSurfaceDataContainer
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

    class ConstantDataContainer
    {
    public:
        ///@name Life cycle
        ///@{

        ConstantDataContainer(
            const Element& rElement,
            const GeometryData::IntegrationMethod& rIntegrationMethod,
            const ProcessInfo& rProcessInfo)
            : mrGeometry(rElement.GetGeometry()),
              mrIntegrationMethod(rIntegrationMethod)
        {
            mHelmholtzRadius = rProcessInfo[HELMHOLTZ_RADIUS];
            const auto& integration_points = mrGeometry.IntegrationPoints(mrIntegrationMethod);

            array_1d<double, 3> normal = ZeroVector(3);

            for (IndexType point_number = 0;
                point_number < integration_points.size(); ++point_number) {
                normal += mrGeometry.UnitNormal(point_number, mrIntegrationMethod);
            }

            normal /= integration_points.size();
            normal /= MathUtils<double>::Norm3(normal);

            const BoundedMatrix<double, 3, 3>& id_matrix = IdentityMatrix(3, 3);

            noalias(mTangentProjectionMatrix) = id_matrix - outer_prod(normal, normal);
        }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        const GeometryType& mrGeometry;

        const GeometryData::IntegrationMethod& mrIntegrationMethod;

        BoundedMatrix<double, 3, 3> mTangentProjectionMatrix;

        double mHelmholtzRadius;

        ///@}
        ///@name Friend classes
        ///@{

        friend class HelmholtzSurfaceDataContainer;

        ///@}
    };

    ///@}
    ///@name Life cycle
    ///@{

    HelmholtzSurfaceDataContainer(const Geometry<Node>& rGeometry)
        : mpSolidGeometry(EntityCalculationUtils::CreateSolidGeometry(rGeometry))
    {
    }

    void AddStiffnessGaussPointContributions(
        Matrix& rStiffnessMatrix,
        const double W,
        const IndexType IntegrationPoint,
        const ConstantDataContainer& rConstantData) const
    {
        Matrix DN_DX;
        EntityCalculationUtils::CalculateSurfaceElementShapeDerivatives(DN_DX, *mpSolidGeometry, rConstantData.mrGeometry, rConstantData.mrIntegrationMethod, IntegrationPoint);
        const Matrix& rdNdX = prod(DN_DX, rConstantData.mTangentProjectionMatrix);

        const BoundedMatrix<double, NumberOfNodes, NumberOfNodes>& A_dirc = W * rConstantData.mHelmholtzRadius * rConstantData.mHelmholtzRadius * prod(rdNdX, trans(rdNdX));

        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const IndexType index_i = i * NumberOfVariables;
            for (IndexType k = 0; k < NumberOfNodes; ++k) {
                const IndexType index_k = k * NumberOfVariables;
                for (IndexType j = 0; j < NumberOfVariables; ++j) {
                    rStiffnessMatrix(index_i + j, index_k + j) += A_dirc(i, k);
                }
            }
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Geometry<Node>::Pointer mpSolidGeometry;

    ///@}
};

} // namespace Kratos