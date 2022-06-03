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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "incompressible_potential_flow_velocity_inlet_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
int IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = BaseType::Check(rCurrentProcessInfo);

    const auto& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_POTENTIAL, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsInlet(*this)) {
        const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);
        KRATOS_ERROR_IF(norm_2(r_normal) == 0.0)
            << "NORMAL is not initialized for inlet condition " << this->Info();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = Condition::GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& ConditionDofList,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (ConditionDofList.size() != TNumNodes) {
        ConditionDofList.resize(TNumNodes);
    }

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        ConditionDofList[i] = Condition::GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::GetValuesVector(
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
            r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::CalculateLocalSystem(
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
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::CalculateLeftHandSide(
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
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    if (IsInlet(*this)) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        CalculateConditionGeometryData(
            this->GetGeometry(), this->GetIntegrationMethod(), gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        array_1d<double, 3> r_normal = this->GetValue(NORMAL);
        const double normal_magnitude = norm_2(r_normal);
        KRATOS_ERROR_IF(normal_magnitude == 0)
            << "NORMAL is not properly initialized for condition id "
            << this->Id() << ".";
        r_normal /= normal_magnitude;

        array_1d<double, 3> velocity;

        for (IndexType g = 0; g < num_gauss_points; ++g) {
            const Vector gauss_shape_functions = row(shape_functions, g);

            FluidCalculationUtilities::EvaluateInPoint(
                this->GetGeometry(), gauss_shape_functions,
                std::tie(velocity, VELOCITY));

            const double velocity_potential_flux =
                inner_prod(velocity, r_normal) * gauss_weights[g];

            noalias(rRightHandSideVector) += gauss_shape_functions * velocity_potential_flux;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePotentialFlowVelocityInletCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "IncompressiblePotentialFlowVelocityInletCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::PrintData(
    std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::save(
    Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>::load(
    Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class IncompressiblePotentialFlowVelocityInletCondition<2, 2>;
template class IncompressiblePotentialFlowVelocityInletCondition<3, 3>;

} // namespace Kratos.
