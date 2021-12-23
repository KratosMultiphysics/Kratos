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

#include "adjoint_diffusion_element.h"
#include "laplacian_element.h"

#include "convection_diffusion_application_variables.h"

#include "includes/checks.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

template<class PrimalElement>
AdjointDiffusionElement<PrimalElement>::AdjointDiffusionElement(IndexType NewId, typename GeometryType::Pointer pGeometry):
    PrimalElement(NewId, pGeometry)
{}

template<class PrimalElement>
AdjointDiffusionElement<PrimalElement>::AdjointDiffusionElement(
    IndexType NewId, typename GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    PrimalElement(NewId, pGeometry, pProperties)
{}

template<class PrimalElement>
AdjointDiffusionElement<PrimalElement>::~AdjointDiffusionElement() {}

template<class PrimalElement>
Element::Pointer AdjointDiffusionElement<PrimalElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointDiffusionElement<PrimalElement>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template<class PrimalElement>
Element::Pointer AdjointDiffusionElement<PrimalElement>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointDiffusionElement<PrimalElement>>(NewId, pGeometry, pProperties);
}

template<class PrimalElement>
void AdjointDiffusionElement<PrimalElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Delegating LHS matrix to base class (the adjoint system matrix is also a laplacian)
    PrimalElement::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Setting the RHS vector to zero
    noalias(rRightHandSideVector) = ZeroVector(rLeftHandSideMatrix.size2());
}

template<class PrimalElement>
void AdjointDiffusionElement<PrimalElement>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const Geometry<Node<3>>& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rRightHandSideVector.size() != num_nodes)
    {
        rRightHandSideVector.resize(num_nodes,false);
    }

    noalias(rRightHandSideVector) = ZeroVector(num_nodes);
}

template<class PrimalElement>
void AdjointDiffusionElement<PrimalElement>::GetValuesVector(Vector& rValues, int Step) const 
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
void AdjointDiffusionElement<PrimalElement>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
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
void AdjointDiffusionElement<PrimalElement>::GetDofList(
    DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rElementalDofList.size() != num_nodes)
    {
        rElementalDofList.resize(num_nodes);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rElementalDofList[i] = r_geom[i].pGetDof(ADJOINT_HEAT_TRANSFER);
    }
}

template<class PrimalElement>
int AdjointDiffusionElement<PrimalElement>::Check(const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    // We check for the primal element too
    // (I am not calling PrimalElement::Check directly because I do not need the Primal Unknown to be a DOF)
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    const auto& r_geom = this->GetGeometry();

    for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
    {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_HEAT_TRANSFER, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_HEAT_TRANSFER, r_node);

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusivity_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_volume_source_var, r_node);
    }

    KRATOS_CATCH("")
    return Element::Check(rProcessInfo);
}

template<class PrimalElement>
std::string AdjointDiffusionElement<PrimalElement>::Info() const
{
    std::stringstream buffer;
    buffer << "AdjointDiffusionElement #" << this->Id();
    return buffer.str();
}

template<class PrimalElement>
void AdjointDiffusionElement<PrimalElement>::PrintInfo(std::ostream& rOStream) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int num_nodes = r_geom.PointsNumber();
    rOStream << "AdjointDiffusionElement" << dimension << "D" << num_nodes << "N";
}

template<class PrimalElement>
void AdjointDiffusionElement<PrimalElement>::CalculateSensitivityMatrix(
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
    Vector diffusivities = ZeroVector(num_nodes);
    Vector volume_source = ZeroVector(num_nodes);

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_primal_values_variable = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_variable = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_variable = r_settings.GetVolumeSourceVariable();

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        primal_values[i] = r_geom[i].FastGetSolutionStepValue(r_primal_values_variable);
        diffusivities[i] = r_geom[i].FastGetSolutionStepValue(r_diffusivity_variable);
        volume_source[i] = r_geom[i].FastGetSolutionStepValue(r_volume_source_variable);
    }

    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        Matrix shape_function_local_gradients(num_nodes,dimension);
        Matrix shape_function_global_gradients(num_nodes,dimension);
        Matrix jacobian(dimension,dimension);
        Matrix jacobian_inv(dimension,dimension);

        const Matrix& N_values = r_geom.ShapeFunctionsValues(integration_method);

        for (unsigned int g = 0; g < num_integration_points; g++)
        {
            noalias(shape_function_local_gradients) = r_geom.ShapeFunctionLocalGradient(g, integration_method);
            r_geom.Jacobian(jacobian, g, integration_method);
            GeometricalSensitivityUtility geometrical_sensitivity_utility(jacobian,shape_function_local_gradients);

            double det_j;
            MathUtils<double>::GeneralizedInvertMatrix(jacobian, jacobian_inv, det_j);
            noalias(shape_function_global_gradients) = prod(shape_function_local_gradients, jacobian_inv);
            const double weight = integration_points[g].Weight();
            const auto N = row(N_values, g);
            const double diffusivity = inner_prod(N, diffusivities);
            const double volume_source_gauss = inner_prod(N, volume_source);

            Vector primal_gradient = prod(trans(shape_function_global_gradients),primal_values);
            Matrix laplacian_operator = prod(shape_function_global_gradients, trans(shape_function_global_gradients));
            Vector laplacian_rhs = prod(laplacian_operator, primal_values);

            for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                double det_j_deriv;
                GeometricalSensitivityUtility::ShapeFunctionsGradientType shape_function_gradient_deriv;
                geometrical_sensitivity_utility.CalculateSensitivity(deriv, det_j_deriv, shape_function_gradient_deriv);

                Vector aux = prod(trans(shape_function_gradient_deriv), primal_values);

                // Calculation by product rule of d/dX_l ( w * J * dN_i/dx_k * dN_j/dx_k * Temperature_j )
                // Plus source term contribution d/dX_l ( w * J * N_i * N_j * VolumeSource_j )
                for (unsigned int i = 0; i < num_nodes; i++)
                {
                    // partial derivative of the determinant of Jacobian: w * conductivity * dJ/dX * dN_i/dx_k * dN_j/dx_k * Temperature_j
                    double contribution_li = det_j_deriv * laplacian_rhs[i];

                    for (unsigned int k = 0; k < dimension; k++)
                    {
                        // partial derivative of the "test function" shape function gradient w * J * d/dX(dN_i/dx_k) * dN_j/dx_k * Temperature_j
                        contribution_li += det_j * shape_function_gradient_deriv(i,k) * primal_gradient[k];

                        // partial derivative of the "primal variable" shape function gradient w * J * dN_i/dx_k * d/dX(dN_j/dx_k) * Temperature_j
                        contribution_li += det_j * shape_function_global_gradients(i,k) * aux[k];
                    }

                    contribution_li *= weight*diffusivity;
                    // plus partial derivative of source term: w * dJ/DX * N_i * N_j * VolumeSource_j
                    contribution_li -= weight * det_j_deriv * N[i] * volume_source_gauss;

                    rOutput(deriv.NodeIndex * dimension + deriv.Direction, i) += contribution_li;
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

template class AdjointDiffusionElement<LaplacianElement>;

}