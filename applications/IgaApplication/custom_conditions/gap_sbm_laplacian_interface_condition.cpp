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
#include "includes/convection_diffusion_settings.h"
#include "custom_conditions/gap_sbm_laplacian_interface_condition.h"

namespace Kratos
{

void GapSbmLaplacianInterfaceCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}

void GapSbmLaplacianInterfaceCondition::InitializeMemberVariables()
{
    const auto& r_geometry = GetGeometry();
    const auto& r_surrogate_plus = GetGeometryPlus();
    const auto& r_DN_De = r_surrogate_plus.ShapeFunctionsLocalGradients(r_surrogate_plus.GetDefaultIntegrationMethod());

    mDim = r_DN_De[0].size2();
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    mBasisFunctionsOrder *= 2;

    // normal from true geometry
    array_1d<double,3> normal_ps = r_geometry.Normal(0, GetIntegrationMethod());
    normal_ps /= MathUtils<double>::Norm(normal_ps);
    mNormalPhysicalSpace = normal_ps;
    SetValue(NORMAL, mNormalPhysicalSpace);

    // integration weight
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const double integration_weight = r_integration_points[0].Weight();
    SetValue(INTEGRATION_WEIGHT, integration_weight);

    // penalty setup
    double penalty = GetProperties()[PENALTY_FACTOR];
    mNitschePenalty = 1.0;
    if (penalty == -1.0) {
        mPenalty = 0.0;
        mNitschePenalty = -1.0; // penalty-free
    } else {
        Vector mesh_size_uv = ZeroVector(3);
        if (Has(KNOT_SPAN_SIZES)) mesh_size_uv = GetValue(KNOT_SPAN_SIZES);
        double h = 1.0;
        if (mesh_size_uv.size() >= 2) {
            h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
            if (mDim == 3 && mesh_size_uv.size() >= 3) h = std::min(h, mesh_size_uv[2]);
        }
        mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * penalty / h;
    }
}

void GapSbmLaplacianInterfaceCondition::InitializeSbmMemberVariables()
{
    const auto& r_true = GetGeometry();
    const auto& r_surrogate_plus = GetGeometryPlus();
    const auto& r_surrogate_minus = GetGeometryMinus();

    mDistanceVectorPlus.resize(3);
    noalias(mDistanceVectorPlus) = r_true.Center().Coordinates() - r_surrogate_plus.Center().Coordinates();
    mDistanceVectorMinus.resize(3);
    noalias(mDistanceVectorMinus) = r_true.Center().Coordinates() - r_surrogate_minus.Center().Coordinates();
}

void GapSbmLaplacianInterfaceCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();
    const std::size_t number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;
    const std::size_t mat_size = number_of_control_points * 1;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
    
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void GapSbmLaplacianInterfaceCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();

    const auto& r_plus = GetGeometryPlus();
    const auto& r_minus = GetGeometryMinus();

    const SizeType n_plus = r_plus.PointsNumber();
    const SizeType n_minus = r_minus.PointsNumber();
    const SizeType n_total = n_plus + n_minus;

    if (rLeftHandSideMatrix.size1() != n_total || rLeftHandSideMatrix.size2() != n_total) {
        rLeftHandSideMatrix.resize(n_total, n_total, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_total, n_total);

    const double weight = GetValue(INTEGRATION_WEIGHT);
    const double k_plus = GetProperties().GetValue(r_diffusivity_var); // same props on both sides by default
    const double k_minus = GetProperties().GetValue(r_diffusivity_var);

    // N_sum
    Vector N_sum_plus(n_plus);
    ComputeTaylorExpansionContribution(r_plus, mDistanceVectorPlus, N_sum_plus);
    Vector N_sum_minus(n_minus);
    ComputeTaylorExpansionContribution(r_minus, mDistanceVectorMinus, N_sum_minus);

    // grad(H_sum)
    Matrix grad_H_plus_T(3, n_plus);
    ComputeGradientTaylorExpansionContribution(r_plus, mDistanceVectorPlus, grad_H_plus_T);
    Matrix grad_H_plus = trans(grad_H_plus_T); // [n_plus x 3]

    Matrix grad_H_minus_T(3, n_minus);
    ComputeGradientTaylorExpansionContribution(r_minus, mDistanceVectorMinus, grad_H_minus_T);
    Matrix grad_H_minus = trans(grad_H_minus_T); // [n_minus x 3]

    // DN dot n on each side
    Vector DNn_plus(n_plus, 0.0);
    for (IndexType a = 0; a < n_plus; ++a)
        for (IndexType d = 0; d < mDim; ++d)
            DNn_plus[a] += grad_H_plus(a, d) * mNormalPhysicalSpace[d];
    Vector DNn_minus(n_minus, 0.0);
    for (IndexType b = 0; b < n_minus; ++b)
        for (IndexType d = 0; d < mDim; ++d)
            DNn_minus[b] += grad_H_minus(b, d) * mNormalPhysicalSpace[d];

    const IndexType off = n_plus;

    // Consistency terms: -[[v]] * {{∂n u}}
    for (IndexType i = 0; i < n_plus; ++i) {
        for (IndexType j = 0; j < n_plus; ++j)
            rLeftHandSideMatrix(i, j) -= 0.5 * weight * N_sum_plus(i) * k_plus * DNn_plus[j];
        for (IndexType j = 0; j < n_minus; ++j)
            rLeftHandSideMatrix(i, off + j) -= 0.5 * weight * N_sum_plus(i) * k_minus * DNn_minus[j];
    }
    // rows: minus side + sign
    for (IndexType i = 0; i < n_minus; ++i) {
        for (IndexType j = 0; j < n_plus; ++j)
            rLeftHandSideMatrix(off + i, j) += 0.5 * weight * N_sum_minus(i) * k_plus * DNn_plus[j];
        for (IndexType j = 0; j < n_minus; ++j)
            rLeftHandSideMatrix(off + i, off + j) += 0.5 * weight * N_sum_minus(i) * k_minus * DNn_minus[j];
    }

    // Symmetric Nitsche terms: -{{∂n v}} * [[u]]
    for (IndexType i = 0; i < n_plus; ++i) {
        for (IndexType j = 0; j < n_plus; ++j)
            rLeftHandSideMatrix(i, j) -= mNitschePenalty * 0.5 * weight * k_plus * DNn_plus[i] * N_sum_plus(j);
        for (IndexType j = 0; j < n_minus; ++j)
            rLeftHandSideMatrix(i, off + j) += mNitschePenalty * 0.5 * weight * k_minus * DNn_plus[i] * N_sum_minus(j);
    }
    for (IndexType i = 0; i < n_minus; ++i) {
        for (IndexType j = 0; j < n_plus; ++j)
            rLeftHandSideMatrix(off + i, j) -= mNitschePenalty * 0.5 * weight * k_plus * DNn_minus[i] * N_sum_plus(j);
        for (IndexType j = 0; j < n_minus; ++j)
            rLeftHandSideMatrix(off + i, off + j) += mNitschePenalty * 0.5 * weight * k_minus * DNn_minus[i] * N_sum_minus(j);
    }

    // Penalty on jump: γ [[v]] * [[u]]
    if (mPenalty != 0.0) {
        for (IndexType i = 0; i < n_plus; ++i) {
            for (IndexType j = 0; j < n_plus; ++j)
                rLeftHandSideMatrix(i, j) += weight * mPenalty * N_sum_plus(i) * N_sum_plus(j);
            for (IndexType j = 0; j < n_minus; ++j)
                rLeftHandSideMatrix(i, off + j) -= weight * mPenalty * N_sum_plus(i) * N_sum_minus(j);
        }
        for (IndexType i = 0; i < n_minus; ++i) {
            for (IndexType j = 0; j < n_plus; ++j)
                rLeftHandSideMatrix(off + i, j) -= weight * mPenalty * N_sum_minus(i) * N_sum_plus(j);
            for (IndexType j = 0; j < n_minus; ++j)
                rLeftHandSideMatrix(off + i, off + j) += weight * mPenalty * N_sum_minus(i) * N_sum_minus(j);
        }
    }
}

void GapSbmLaplacianInterfaceCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType lhs;
    CalculateLeftHandSide(lhs, rCurrentProcessInfo);

    const SizeType size = lhs.size1();
    if (rRightHandSideVector.size() != size)
        rRightHandSideVector.resize(size, false);

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_plus = GetGeometryPlus();
    const auto& r_minus = GetGeometryMinus();
    Vector values(r_plus.PointsNumber() + r_minus.PointsNumber());
    for (IndexType i = 0; i < r_plus.PointsNumber(); ++i) values[i] = r_plus[i].FastGetSolutionStepValue(r_unknown_var);
    for (IndexType i = 0; i < r_minus.PointsNumber(); ++i) values[r_plus.PointsNumber() + i] = r_minus[i].FastGetSolutionStepValue(r_unknown_var);

    noalias(rRightHandSideVector) = -prod(lhs, values);
}

void GapSbmLaplacianInterfaceCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_plus = GetGeometryPlus();
    const auto& r_minus = GetGeometryMinus();
    const SizeType n_plus = r_plus.PointsNumber();
    const SizeType n_minus = r_minus.PointsNumber();
    const SizeType n_total = n_plus + n_minus;
    if (rResult.size() != n_total) rResult.resize(n_total, false);
    for (IndexType i = 0; i < n_plus; ++i) rResult[i] = r_plus[i].GetDof(r_unknown_var).EquationId();
    for (IndexType i = 0; i < n_minus; ++i) rResult[n_plus + i] = r_minus[i].GetDof(r_unknown_var).EquationId();
}

void GapSbmLaplacianInterfaceCondition::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_plus = GetGeometryPlus();
    const auto& r_minus = GetGeometryMinus();
    const SizeType n_plus = r_plus.PointsNumber();
    const SizeType n_minus = r_minus.PointsNumber();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(n_plus + n_minus);
    for (IndexType i = 0; i < n_plus; ++i) rElementalDofList.push_back(r_plus[i].pGetDof(r_unknown_var));
    for (IndexType i = 0; i < n_minus; ++i) rElementalDofList.push_back(r_minus[i].pGetDof(r_unknown_var));
}

void GapSbmLaplacianInterfaceCondition::ComputeTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Vector& H_sum_vec)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    if (H_sum_vec.size() != number_of_control_points) H_sum_vec.resize(number_of_control_points);
    // derivatives up to order mBasisFunctionsOrder
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term = 0.0;
        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    const double derivative = r_shape_function_derivatives(i, k);
                    H_taylor_term += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
                int countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        const double derivative = r_shape_function_derivatives(i, countDerivativeId);
                        H_taylor_term += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0, i);
    }
}

void GapSbmLaplacianInterfaceCondition::ComputeGradientTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Matrix& grad_H_sum)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(rGeometry.GetDefaultIntegrationMethod());
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points) grad_H_sum.resize(3, number_of_control_points);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term_X = 0.0, H_taylor_term_Y = 0.0, H_taylor_term_Z = 0.0;
        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n - 1];
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k);
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k + 1);
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n - 1];
                IndexType countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        const double derivative = shapeFunctionDerivatives(i, countDerivativeId);
                        if (k_x >= 1) H_taylor_term_X += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x - 1, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        if (k_y >= 1) H_taylor_term_Y += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y - 1, rDistanceVector[2], k_z);
                        if (k_z >= 1) H_taylor_term_Z += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z - 1);
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

double GapSbmLaplacianInterfaceCondition::ComputeTaylorTerm(const double derivative, const double dx, const IndexType n_k, const double dy, const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

double GapSbmLaplacianInterfaceCondition::ComputeTaylorTerm3D(const double derivative, const double dx, const IndexType k_x, const double dy, const IndexType k_y, const double dz, const IndexType k_z)
{
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
}

} // namespace Kratos

