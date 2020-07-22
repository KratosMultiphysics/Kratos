//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "incompressible_potential_flow_pressure_body_force_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
int IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::Check(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = BaseType::Check(rCurrentProcessInfo);

    const auto& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_POTENTIAL, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_node);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE_POTENTIAL, r_node);
    }

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY;

    const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);
    KRATOS_ERROR_IF(norm_2(r_normal) == 0.0)
        << "NORMAL is not initialized for " << this->Info();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = Condition::GetGeometry()[i].GetDof(PRESSURE_POTENTIAL).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& ConditionDofList,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (ConditionDofList.size() != TNumNodes) {
        ConditionDofList.resize(TNumNodes);
    }

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        ConditionDofList[i] = Condition::GetGeometry()[i].pGetDof(PRESSURE_POTENTIAL);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TNumNodes) {
        rValues.resize(TNumNodes, false);
    }

    const auto& r_geometry = this->GetGeometry();
    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rValues[local_index++] =
            r_geometry[i_node].FastGetSolutionStepValue(PRESSURE_POTENTIAL, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Calculate LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    const auto& r_geometry = this->GetGeometry();
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    RansCalculationUtilities::CalculateConditionGeometryData(
        r_geometry, this->GetIntegrationMethod(), gauss_weights, shape_functions);
    const IndexType num_gauss_points = gauss_weights.size();

    array_1d<double, 3> r_normal = this->GetValue(NORMAL);
    r_normal = r_normal * (1.0 / norm_2(r_normal));

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Vector gauss_shape_functions = row(shape_functions, g);

        const array_1d<double, 3>& r_body_force = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, BODY_FORCE, gauss_shape_functions);
        const double density = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, DENSITY, gauss_shape_functions);

        const double value =
            density * inner_prod(r_body_force, r_normal) * gauss_weights[g];

        for (IndexType a = 0; a < TNumNodes; ++a) {
            rRightHandSideVector[a] += value * gauss_shape_functions[a];
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePotentialFlowPressureBodyForceCondition" << TDim << "D"
           << " ID: " << this->Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "IncompressiblePotentialFlowPressureBodyForceCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::PrintData(
    std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureBodyForceCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class IncompressiblePotentialFlowPressureBodyForceCondition<2, 2>;
template class IncompressiblePotentialFlowPressureBodyForceCondition<3, 3>;

} // namespace Kratos.
