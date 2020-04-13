//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
#include "incompressible_potential_flow_velocity_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>& IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::operator=(
    IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityCondition>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityCondition>(
        NewId, pGeom, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition = Create(
        NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
int IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = BaseType::Check(rCurrentProcessInfo);

    const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);
    KRATOS_ERROR_IF(norm_2(r_normal) == 0.0)
        << "NORMAL is not initialized for " << this->Info();

    const GeometryType& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_POTENTIAL, r_node);
    }

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; ++i)
        rResult[i] = Condition::GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo)
{
    if (ConditionDofList.size() != TNumNodes)
        ConditionDofList.resize(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; ++i)
        ConditionDofList[i] = Condition::GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::GetValuesVector(
    VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] =
            rGeom[iNode].FastGetSolutionStepValue(VELOCITY_POTENTIAL, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Calculate LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes, false);

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    if (RansCalculationUtilities::IsInlet(*this))
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            this->GetGeometry(), this->GetIntegrationMethod(), gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        array_1d<double, 3> r_normal = this->GetValue(NORMAL);
        r_normal = r_normal * (1.0 / norm_2(r_normal));

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Vector gauss_shape_functions = row(shape_functions, g);

            const array_1d<double, 3>& r_velocity = RansCalculationUtilities::EvaluateInPoint(
                this->GetGeometry(), VELOCITY, gauss_shape_functions);

            const double velocity_potential_flux =
                inner_prod(r_velocity, r_normal) * gauss_weights[g];

            noalias(rRightHandSideVector) += gauss_shape_functions * velocity_potential_flux;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePotentialFlowVelocityCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "IncompressiblePotentialFlowVelocityCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class IncompressiblePotentialFlowVelocityCondition<2, 2>;
template class IncompressiblePotentialFlowVelocityCondition<3, 3>;

} // namespace Kratos.
