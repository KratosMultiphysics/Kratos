// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#include "adjoint_flux_condition.h"
#include "flux_condition.h"

#include "convection_diffusion_application_variables.h"

#include "includes/checks.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

template<class PrimalCondition>
AdjointFluxCondition<PrimalCondition>::AdjointFluxCondition(IndexType NewId, typename GeometryType::Pointer pGeometry):
    PrimalCondition(NewId, pGeometry)
{}

template<class PrimalCondition>
AdjointFluxCondition<PrimalCondition>::AdjointFluxCondition(
    IndexType NewId, typename GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    PrimalCondition(NewId, pGeometry, pProperties)
{}

template<class PrimalCondition>
AdjointFluxCondition<PrimalCondition>::~AdjointFluxCondition() {}

template<class PrimalCondition>
Condition::Pointer AdjointFluxCondition<PrimalCondition>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointFluxCondition<PrimalCondition>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template<class PrimalCondition>
Condition::Pointer AdjointFluxCondition<PrimalCondition>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointFluxCondition<PrimalCondition>>(NewId, pGeometry, pProperties);
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // Delegating LHS matrix to base class (the adjoint system matrix is also a laplacian)
    PrimalCondition::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Setting the RHS vector to zero
    noalias(rRightHandSideVector) = ZeroVector(rLeftHandSideMatrix.size2());
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const Geometry<Node<3>>& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rRightHandSideVector.size() != num_nodes)
    {
        rRightHandSideVector.resize(num_nodes,false);
    }

    noalias(rRightHandSideVector) = ZeroVector(num_nodes);
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::GetValuesVector(Vector& rValues, int Step)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rValues.size() != num_nodes)
    {
        rValues.resize(num_nodes,false);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rValues[i] = r_geom[i].FastGetSolutionStepValue(ADJOINT_HEAT_TRANSFER, Step);
    }
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rResult.size() != num_nodes)
    {
        rResult.resize(num_nodes,false);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rResult[i] = r_geom[i].GetDof(ADJOINT_HEAT_TRANSFER).EquationId();
    }
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::GetDofList(
    DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rConditionDofList.size() != num_nodes)
    {
        rConditionDofList.resize(num_nodes);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rConditionDofList[i] = r_geom[i].pGetDof(ADJOINT_HEAT_TRANSFER);
    }
}

template<class PrimalCondition>
int AdjointFluxCondition<PrimalCondition>::Check(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();
    for (unsigned int i = 0; i < num_nodes; i++)
    {
        const Node<3>& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_HEAT_TRANSFER, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_HEAT_TRANSFER, r_node);
    }

    return PrimalCondition::Check(rProcessInfo);
    KRATOS_CATCH("")
}

template<class PrimalCondition>
std::string AdjointFluxCondition<PrimalCondition>::Info() const
{
    std::stringstream buffer;
    buffer << "AdjointFluxCondition #" << this->Id();
    return buffer.str();
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::PrintInfo(std::ostream& rOStream) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int num_nodes = r_geom.PointsNumber();
    rOStream << "AdjointFluxCondition" << dimension << "D" << num_nodes << "N";
}

template<class PrimalCondition>
void AdjointFluxCondition<PrimalCondition>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int num_nodes = r_geom.PointsNumber();
    const unsigned int sensitivity_size = dimension * num_nodes;

    if (rOutput.size1() != sensitivity_size || rOutput.size2() != num_nodes)
    {
        rOutput.resize(sensitivity_size,num_nodes,false);
    }
    noalias(rOutput) = ZeroMatrix(sensitivity_size,num_nodes);

    const auto integration_method = this->GetIntegrationMethod();
    const auto integration_points = r_geom.IntegrationPoints(integration_method);
    const unsigned int num_integration_points = integration_points.size();
    Vector primal_values = ZeroVector(num_nodes);
    PrimalCondition::GetValuesVector(primal_values);

    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        Matrix shape_function_local_gradients(num_nodes,dimension);
        Matrix shape_function_global_gradients(num_nodes,dimension);
        Matrix jacobian(dimension,dimension);
        Matrix jacobian_inv(dimension,dimension);

        for (unsigned int g = 0; g < num_integration_points; g++)
        {

        }
    }
    else
    {
        KRATOS_ERROR << "Unsupported design variable " << rDesignVariable << std::endl;
    }

    KRATOS_CATCH("")
}

template class AdjointFluxCondition<FluxCondition<2>>;
template class AdjointFluxCondition<FluxCondition<3>>;

}