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
    const Geometry<Node<3>>& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rLeftHandSideMatrix.size1() != num_nodes || rLeftHandSideMatrix.size2() != num_nodes)
    {
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,0);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(num_nodes,num_nodes);

    if (rRightHandSideVector.size() != num_nodes)
    {
        rRightHandSideVector.resize(num_nodes,false);
    }

    noalias(rRightHandSideVector) = ZeroVector(num_nodes);
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
    Vector nodal_flux = ZeroVector(num_nodes);
    for (unsigned int i = 0; i < num_nodes; i++)
    {
        nodal_flux[i] = r_geom[i].FastGetSolutionStepValue(HEAT_FLUX);
    }

    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        Matrix shape_function_local_gradients(num_nodes,dimension);
        Matrix jacobian(dimension,dimension);
        Matrix jacobian_inv(dimension,dimension);

        Matrix shape_functions = r_geom.ShapeFunctionsValues(integration_method);

        for (unsigned int g = 0; g < num_integration_points; g++)
        {
            noalias(shape_function_local_gradients) = r_geom.ShapeFunctionLocalGradient(g, integration_method);
            r_geom.Jacobian(jacobian, g, integration_method);
            GeometricalSensitivityUtility geometrical_sensitivity_utility(jacobian,shape_function_local_gradients);

            Vector N = row(shape_functions, g);
            double q_gauss = prod(N,nodal_flux);

            double det_j;
            MathUtils<double>::GeneralizedInvertMatrix(jacobian, jacobian_inv, det_j);
            const double weight = integration_points[g].Weight();

            Vector primal_gradient = prod(trans(shape_function_global_gradients),primal_values);
            Matrix laplacian_operator = prod(shape_function_global_gradients, trans(shape_function_global_gradients));
            Vector laplacian_rhs = prod(laplacian_operator, primal_values);

            for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                const unsigned int l = deriv.NodeIndex * dimension + deriv.Direction;
                double det_j_deriv;
                GeometricalSensitivityUtility::ShapeFunctionsGradientType shape_function_gradient_deriv;
                geometrical_sensitivity_utility.CalculateSensitivity(deriv, det_j_deriv, shape_function_gradient_deriv);

                Vector aux = prod(trans(shape_function_gradient_deriv), primal_values);

                // d/dX_l (w * J * N_i * N_j * q_j) = w * N_i * N_j * q_j * dJ/dX_l
                // Note that N_j * q_j = q_gauss
                for (unsigned int i = 0; i < num_nodes; i++)
                {
                    rOutput(deriv.NodeIndex * dimension + deriv.Direction, i) += weight * N[i] * q_gauss * det_j_deriv;
                }
            }
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