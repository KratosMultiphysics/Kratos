//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_elements/gap_sbm_laplacian_element.h"

namespace Kratos
{

GapSbmLaplacianElement::GapSbmLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry)
: Element(NewId, pGeometry)
{
}

GapSbmLaplacianElement::GapSbmLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
: Element(NewId, pGeometry, pProperties)
{
}

GapSbmLaplacianElement::~GapSbmLaplacianElement()
{
}

Element::Pointer GapSbmLaplacianElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmLaplacianElement>(NewId, GetSurrogateGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer GapSbmLaplacianElement::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmLaplacianElement>(NewId, pGeom, pProperties);
}

void GapSbmLaplacianElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_WATCH("gap interface laplacian")
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    KRATOS_WATCH("end")

    // Set a single-point integration weight (as in standard Laplacian)
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    SetValue(INTEGRATION_WEIGHT, r_integration_points[0].Weight());
}

void GapSbmLaplacianElement::InitializeMemberVariables()
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_surrogate_geometry.ShapeFunctionsLocalGradients(r_surrogate_geometry.GetDefaultIntegrationMethod());

    mDim = r_DN_De[0].size2();

    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1; // heuristic consistent with other GAP-SBM elements
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    mBasisFunctionsOrder *= 2; // follow existing GAP-SBM elements convention
}

void GapSbmLaplacianElement::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();
}

void GapSbmLaplacianElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const std::size_t mat_size = GetSurrogateGeometry().size() * 1;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void GapSbmLaplacianElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // e.g., conductivity

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const unsigned int number_of_points = r_surrogate_geometry.PointsNumber();

    const double weight = GetValue(INTEGRATION_WEIGHT);
    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);

    if (rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points, number_of_points);

    // compute Taylor expansion contribution: grad(H_sum)
    Matrix grad_H_sum_T(3, number_of_points);
    ComputeGradientTaylorExpansionContribution(grad_H_sum_T);
    Matrix grad_H_sum = trans(grad_H_sum_T); // size: [n_points x 3]

    noalias(rLeftHandSideMatrix) += weight * conductivity * prod(grad_H_sum, trans(grad_H_sum));
}

void GapSbmLaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const unsigned int number_of_points = r_surrogate_geometry.PointsNumber();

    const double weight = GetValue(INTEGRATION_WEIGHT);
    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);

    if (rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points, false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points);

    // N_sum (Taylor-expanded shape function)
    Vector H_sum_vec(number_of_points);
    ComputeTaylorExpansionContribution(H_sum_vec);

    // grad(H_sum)
    Matrix grad_H_sum_T(3, number_of_points);
    ComputeGradientTaylorExpansionContribution(grad_H_sum_T);
    Matrix grad_H_sum = trans(grad_H_sum_T); // [n_points x 3]

    // Unknown at nodes (on surrogate geometry)
    Vector unknown_old(number_of_points);
    for (IndexType i = 0; i < number_of_points; ++i) {
        unknown_old[i] = r_surrogate_geometry[i].GetSolutionStepValue(r_unknown_var);
    }

    const double source = this->GetValue(r_volume_source_var);

    // RHS = -K*u + source * N
    noalias(rRightHandSideVector) -= weight * conductivity * prod(prod(grad_H_sum, trans(grad_H_sum)), unknown_old);
    noalias(rRightHandSideVector) += weight * H_sum_vec * source;
}

void GapSbmLaplacianElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_geometry = GetSurrogateGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
    }
}

void GapSbmLaplacianElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_geometry = GetSurrogateGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    if (rElementalDofList.size() != number_of_nodes)
        rElementalDofList.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        rElementalDofList[i] = r_geometry[i].pGetDof(r_unknown_var);
    }
}

int GapSbmLaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const auto& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    const auto& r_geom = GetSurrogateGeometry();
    for (unsigned int i = 0; i < r_geom.PointsNumber(); ++i) {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusivity_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_volume_source_var, r_node);
        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
    }

    return Element::Check(rCurrentProcessInfo);
}

Element::IntegrationMethod GapSbmLaplacianElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

// Taylor expansion helpers (copied from GAP-SBM elements for consistency)
void GapSbmLaplacianElement::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
    const auto& r_N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (H_sum_vec.size() != number_of_control_points)
        H_sum_vec.resize(number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term = 0.0;
        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    double derivative = r_shape_function_derivatives(i, k);
                    H_taylor_term += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                int countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i, countDerivativeId);
                        H_taylor_term += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0, i);
    }
}

void GapSbmLaplacianElement::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points)
        grad_H_sum.resize(3, number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term_X = 0.0;
        double H_taylor_term_Y = 0.0;
        double H_taylor_term_Z = 0.0;

        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative_x = shapeFunctionDerivatives(i, k);
                    H_taylor_term_X += ComputeTaylorTerm(derivative_x, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative_y = shapeFunctionDerivatives(i, k + 1);
                    H_taylor_term_Y += ComputeTaylorTerm(derivative_y, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                IndexType countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        double derivative = shapeFunctionDerivatives(i, countDerivativeId);
                        if (k_x >= 1) {
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x - 1, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y - 1, mDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z - 1);
                        }
                        countDerivativeId++;
                    }
                }
            }
        }

        grad_H_sum(0, i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1, i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2, i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else
            grad_H_sum(2, i) = 0.0;
    }
}

double GapSbmLaplacianElement::ComputeTaylorTerm(const double derivative, const double dx, const IndexType n_k, const double dy, const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

double GapSbmLaplacianElement::ComputeTaylorTerm3D(const double derivative, const double dx, const IndexType k_x, const double dy, const IndexType k_y, const double dz, const IndexType k_z)
{
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
}

} // namespace Kratos

