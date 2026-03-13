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
#include "custom_conditions/gap_sbm_laplacian_condition.h"

namespace Kratos
{

void GapSbmLaplacianCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}

void GapSbmLaplacianCondition::InitializeMemberVariables()
{
    const auto& r_true_geometry = GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_surrogate_geometry.ShapeFunctionsLocalGradients(r_surrogate_geometry.GetDefaultIntegrationMethod());

    mDim = r_DN_De[0].size2();
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    mBasisFunctionsOrder *= 2;

    // normal from true geometry
    array_1d<double,3> normal_ps = r_true_geometry.Normal(0, GetIntegrationMethod());
    normal_ps /= MathUtils<double>::Norm(normal_ps);
    mNormalPhysicalSpace = normal_ps;
    SetValue(NORMAL, mNormalPhysicalSpace);

    // integration weight
    const auto& r_integration_points = r_true_geometry.IntegrationPoints(this->GetIntegrationMethod());
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

void GapSbmLaplacianCondition::InitializeSbmMemberVariables()
{
    const auto& r_true = GetGeometry();
    const auto& r_surrogate = GetSurrogateGeometry();
    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_true.Center().Coordinates() - r_surrogate.Center().Coordinates();
}

void GapSbmLaplacianCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const SizeType mat_size = r_surrogate_geometry.PointsNumber();
    
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void GapSbmLaplacianCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const SizeType n = r_surrogate_geometry.PointsNumber();
    if (rLeftHandSideMatrix.size1() != n || rLeftHandSideMatrix.size2() != n) rLeftHandSideMatrix.resize(n, n, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n, n);

    const double weight = GetValue(INTEGRATION_WEIGHT);

    // H_sum and grad(H_sum)
    Vector H_sum_vec(n);
    ComputeTaylorExpansionContribution(H_sum_vec);
    Matrix grad_H_sum_T(3, n);
    ComputeGradientTaylorExpansionContribution(grad_H_sum_T);
    Matrix grad_H_sum = trans(grad_H_sum_T);

    // DN dot n
    Vector DNn(n, 0.0);
    for (IndexType i = 0; i < n; ++i)
        for (IndexType d = 0; d < mDim; ++d)
            DNn[i] += grad_H_sum(i, d) * mNormalPhysicalSpace[d];

    // Assemble Nitsche Dirichlet terms using GAP-SBM expansions
    // -(w, ∂n u)  - mNitschePenalty (∂n w, u)  + γ (w, u)
    // Consistency term: -(w, ∂n u) -> outer product of H_sum_vec and DNn
    noalias(rLeftHandSideMatrix) -= weight * outer_prod(H_sum_vec, DNn);
    // Symmetric Nitsche term: (∂n w, u)
    noalias(rLeftHandSideMatrix) -= weight * mNitschePenalty * outer_prod(DNn, H_sum_vec);

    // Penalty term: γ (w, u)
    if (mPenalty != 0.0)
        noalias(rLeftHandSideMatrix) += weight * mPenalty * outer_prod(H_sum_vec, H_sum_vec);
}

void GapSbmLaplacianCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const SizeType n = r_surrogate_geometry.PointsNumber();
    if (rRightHandSideVector.size() != n) rRightHandSideVector.resize(n, false);
    noalias(rRightHandSideVector) = ZeroVector(n);

    const double weight = GetValue(INTEGRATION_WEIGHT);

    // H_sum and grad(H_sum)
    Vector H_sum_vec(n);
    ComputeTaylorExpansionContribution(H_sum_vec);
    Matrix grad_H_sum_T(3, n);
    ComputeGradientTaylorExpansionContribution(grad_H_sum_T);
    Matrix grad_H_sum = trans(grad_H_sum_T);
    Vector DNn(n, 0.0);
    for (IndexType i = 0; i < n; ++i)
        for (IndexType d = 0; d < mDim; ++d)
            DNn[i] += grad_H_sum(i, d) * mNormalPhysicalSpace[d];

    // Boundary value g from condition data for the unknown variable
    const double g_value = this->GetValue(r_unknown_var);

    // RHS terms corresponding to LHS
    // + (w, k ∂n u) at current u
    Vector u_old(n);
    for (IndexType i = 0; i < n; ++i) 
        u_old[i] = r_surrogate_geometry[i].FastGetSolutionStepValue(r_unknown_var);
    const double DNn_dot_u = inner_prod(DNn, u_old);
    const double H_sum_dot_u = inner_prod(H_sum_vec, u_old);

    noalias(rRightHandSideVector) += H_sum_vec * DNn_dot_u * weight; // (w, ∂n u)
    noalias(rRightHandSideVector) += mNitschePenalty * DNn * H_sum_dot_u * weight; // (∂n w, u)
    noalias(rRightHandSideVector) -= H_sum_vec * H_sum_dot_u * mPenalty * weight; // -γ (w, u)

    // Data terms with g
    noalias(rRightHandSideVector) -= mNitschePenalty * DNn * g_value * weight;
    noalias(rRightHandSideVector) += H_sum_vec * g_value * mPenalty * weight;
}

int GapSbmLaplacianCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << "No penalty factor (PENALTY_FACTOR) defined in property of GapSbmLaplacianCondition" << std::endl;
    return 0;
}

void GapSbmLaplacianCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType n = r_geometry.PointsNumber();
    if (rResult.size() != n) rResult.resize(n, false);
    for (IndexType i = 0; i < n; ++i) rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
}

void GapSbmLaplacianCondition::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType n = r_geometry.PointsNumber();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(n);
    for (IndexType i = 0; i < n; ++i) rElementalDofList.push_back(r_geometry[i].pGetDof(r_unknown_var));
}

void GapSbmLaplacianCondition::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType n = r_geometry.PointsNumber();
    const Matrix& N = r_geometry.ShapeFunctionsValues();
    if (H_sum_vec.size() != n) H_sum_vec.resize(n);
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType ord = 1; ord <= mBasisFunctionsOrder; ++ord) {
        shape_function_derivatives[ord - 1] = r_geometry.ShapeFunctionDerivatives(ord, 0, this->GetIntegrationMethod());
    }
    for (IndexType i = 0; i < n; ++i) {
        double H_taylor = 0.0;
        if (mDim == 2) {
            for (IndexType ord = 1; ord <= mBasisFunctionsOrder; ++ord) {
                Matrix& deriv = shape_function_derivatives[ord - 1];
                for (IndexType k = 0; k <= ord; ++k) {
                    IndexType n_k = ord - k;
                    const double d = deriv(i, k);
                    H_taylor += ComputeTaylorTerm(d, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType ord = 1; ord <= mBasisFunctionsOrder; ++ord) {
                Matrix& deriv = shape_function_derivatives[ord - 1];
                int count = 0;
                for (IndexType kx = ord; kx >= 0; --kx) {
                    for (IndexType ky = ord - kx; ky >= 0; --ky) {
                        IndexType kz = ord - kx - ky;
                        const double d = deriv(i, count++);
                        H_taylor += ComputeTaylorTerm3D(d, mDistanceVector[0], kx, mDistanceVector[1], ky, mDistanceVector[2], kz);
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor + N(0, i);
    }
}

void GapSbmLaplacianCondition::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType n = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType ord = 1; ord <= mBasisFunctionsOrder; ++ord) {
        shape_function_derivatives[ord - 1] = r_geometry.ShapeFunctionDerivatives(ord, 0, this->GetIntegrationMethod());
    }
    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != n) grad_H_sum.resize(3, n);
    for (IndexType i = 0; i < n; ++i) {
        double tx = 0.0, ty = 0.0, tz = 0.0;
        if (mDim == 2) {
            for (IndexType ord = 2; ord <= mBasisFunctionsOrder; ++ord) {
                Matrix& deriv = shape_function_derivatives[ord - 1];
                for (IndexType k = 0; k <= ord - 1; ++k) {
                    IndexType n_k = ord - 1 - k;
                    const double d = deriv(i, k);
                    tx += ComputeTaylorTerm(d, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= ord - 1; ++k) {
                    IndexType n_k = ord - 1 - k;
                    const double d = deriv(i, k + 1);
                    ty += ComputeTaylorTerm(d, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType ord = 2; ord <= mBasisFunctionsOrder; ++ord) {
                Matrix& deriv = shape_function_derivatives[ord - 1];
                int count = 0;
                for (IndexType kx = ord; kx >= 0; --kx) {
                    for (IndexType ky = ord - kx; ky >= 0; --ky) {
                        IndexType kz = ord - kx - ky;
                        const double d = deriv(i, count++);
                        if (kx >= 1) tx += ComputeTaylorTerm3D(d, mDistanceVector[0], kx - 1, mDistanceVector[1], ky, mDistanceVector[2], kz);
                        if (ky >= 1) ty += ComputeTaylorTerm3D(d, mDistanceVector[0], kx, mDistanceVector[1], ky - 1, mDistanceVector[2], kz);
                        if (kz >= 1) tz += ComputeTaylorTerm3D(d, mDistanceVector[0], kx, mDistanceVector[1], ky, mDistanceVector[2], kz - 1);
                    }
                }
            }
        }
        grad_H_sum(0, i) = tx + r_DN_De[0](i, 0);
        grad_H_sum(1, i) = ty + r_DN_De[0](i, 1);
        grad_H_sum(2, i) = (mDim == 3) ? (tz + r_DN_De[0](i, 2)) : 0.0;
    }
}

double GapSbmLaplacianCondition::ComputeTaylorTerm(const double derivative, const double dx, const IndexType n_k, const double dy, const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

double GapSbmLaplacianCondition::ComputeTaylorTerm3D(const double derivative, const double dx, const IndexType k_x, const double dy, const IndexType k_y, const double dz, const IndexType k_z)
{
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
}

} // namespace Kratos
