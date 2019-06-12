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

#include "adjoint_heat_diffusion_element.h"
#include "laplacian_element.h"

#include "convection_diffusion_application_variables.h"

#include "includes/checks.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

template<class PrimalElement>
AdjointHeatDiffusionElement<PrimalElement>::AdjointHeatDiffusionElement(
    IndexType NewId, typename GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    PrimalElement(NewId, pGeometry, pProperties)
{}

template<class PrimalElement>
AdjointHeatDiffusionElement<PrimalElement>::~AdjointHeatDiffusionElement() {}

template<class PrimalElement>
Element::Pointer AdjointHeatDiffusionElement<PrimalElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointHeatDiffusionElement<PrimalElement>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template<class PrimalElement>
Element::Pointer AdjointHeatDiffusionElement<PrimalElement>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointHeatDiffusionElement<PrimalElement>>(NewId, pGeometry, pProperties);
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const Geometry<Node<3>>& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rLeftHandSideMatrix.size1() != num_nodes || rLeftHandSideMatrix.size2() != num_nodes)
    {
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(num_nodes,num_nodes);
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::GetValuesVector(Vector& rValues, int Step)
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

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::EquationIdVector(
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

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::GetDofList(
    DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rElementalDofList.size() != num_nodes)
    {
        rElementalDofList.resize(num_nodes);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_HEAT_TRANSFER));
    }
}

template<class PrimalElement>
int AdjointHeatDiffusionElement<PrimalElement>::Check(const ProcessInfo& rProcessInfo)
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

    return PrimalElement::Check(rProcessInfo);
    KRATOS_CATCH("")
}

template<class PrimalElement>
std::string AdjointHeatDiffusionElement<PrimalElement>::Info() const
{
    std::stringstream buffer;
    buffer << "AdjointHeatDiffusionElement #" << this->Id();
    return buffer.str();
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::PrintInfo(std::ostream& rOStream) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int num_nodes = r_geom.PointsNumber();
    rOStream << "AdjointHeatDiffusionElement" << dimension << "D" << num_nodes << "N";
}

template<class PrimalElement>
void AdjointHeatDiffusionElement<PrimalElement>::CalculateSensitivityMatrix(
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
    PrimalElement::GetValuesVector(primal_values);

    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        Matrix shape_function_local_gradients(num_nodes,dimension);
        Matrix shape_function_global_gradients(num_nodes,dimension);
        Matrix jacobian(num_nodes,num_nodes);
        Matrix jacobian_inv(num_nodes,num_nodes);

        for (unsigned int g = 0; g < num_integration_points; g++)
        {
            noalias(shape_function_local_gradients) = r_geom.ShapeFunctionLocalGradient(g, integration_method);
            r_geom.Jacobian(jacobian, g, integration_method);
            GeometricalSensitivityUtility geometrical_sensitivity_utility(jacobian,shape_function_local_gradients);

            double det_j;
            MathUtils<double>::GeneralizedInvertMatrix(jacobian, jacobian_inv, det_j);
            noalias(shape_function_global_gradients) = prod(shape_function_local_gradients, jacobian_inv);
            const double weight = integration_points[g].Weight();

            Vector primal_gradient = prod(trans(shape_function_global_gradients),primal_values);

            for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                double det_j_deriv;
                GeometricalSensitivityUtility::ShapeFunctionsGradientType shape_function_gradient_deriv;
                geometrical_sensitivity_utility.CalculateSensitivity(deriv, det_j_deriv, shape_function_gradient_deriv);

                for (unsigned int j = 0; j < num_nodes; j++)
                {
                    rOutput(deriv.NodeIndex * dimension + deriv.Direction, j) += 0.0; // write to the matrix!
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

template class AdjointHeatDiffusionElement<LaplacianElement>;

}