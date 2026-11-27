//  KRATOS Multiphysics - IgaApplication
//  License: BSD License

#include <cmath>
#include "custom_conditions/sbm_solid_green_lagrange_condition.h"

namespace Kratos
{

Condition::Pointer SbmSolidGreenLagrangeCondition::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmSolidGreenLagrangeCondition>(NewId, pGeometry, pProperties);
}

Condition::Pointer SbmSolidGreenLagrangeCondition::Create(
    IndexType NewId, NodesArrayType const& rNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmSolidGreenLagrangeCondition>(
        NewId, GetGeometry().Create(rNodes), pProperties);
}

void SbmSolidGreenLagrangeCondition::CalculateKinematics(
    Matrix& rDN_DX, Matrix& rF, Matrix& rB) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType nodes = r_geometry.size();
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "Unsupported dimension " << mDim << std::endl;
    KRATOS_ERROR_IF(strain_size != (mDim == 2 ? 3 : 6)) << "Incompatible strain size" << std::endl;
    const auto& r_local_gradients = r_geometry.ShapeFunctionsLocalGradients(GetIntegrationMethod());
    rDN_DX.resize(nodes, mDim, false);
    Matrix jacobian;
    CalculateInitialJacobian(r_geometry, jacobian);
    Matrix inverse_jacobian(mDim, mDim);
    double determinant = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, determinant);
    noalias(rDN_DX) = prod(r_local_gradients[0], inverse_jacobian);
    rF = IdentityMatrix(mDim);
    for (IndexType a = 0; a < nodes; ++a) {
        const auto& r_u = r_geometry[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                rF(i, j) += r_u[i] * rDN_DX(a, j);
            }
        }
    }
    const double det_f = MathUtils<double>::Det(rF);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_f) && det_f > 0.0)
        << "Invalid deformation gradient determinant " << det_f << " in condition " << Id() << std::endl;

    rB = ZeroMatrix(strain_size, nodes * mDim);
    for (IndexType a = 0; a < nodes; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType col = a * mDim + i;
            rB(0, col) = rF(i, 0) * rDN_DX(a, 0);
            rB(1, col) = rF(i, 1) * rDN_DX(a, 1);
            if (mDim == 2) {
                rB(2, col) = rF(i, 0) * rDN_DX(a, 1) + rF(i, 1) * rDN_DX(a, 0);
            } else {
                rB(2, col) = rF(i, 2) * rDN_DX(a, 2);
                rB(3, col) = rF(i, 0) * rDN_DX(a, 1) + rF(i, 1) * rDN_DX(a, 0);
                rB(4, col) = rF(i, 1) * rDN_DX(a, 2) + rF(i, 2) * rDN_DX(a, 1);
                rB(5, col) = rF(i, 0) * rDN_DX(a, 2) + rF(i, 2) * rDN_DX(a, 0);
            }
        }
    }
}

void SbmSolidGreenLagrangeCondition::SetConstitutiveParameters(
    const Matrix& rF, double DeterminantF, ConstitutiveLaw::Parameters& rValues,
    ConstitutiveVariables& rVariables, bool ComputeTangent) const
{
    auto& r_options = rValues.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, ComputeTangent);
    rValues.SetDeformationGradientF(rF);
    rValues.SetDeterminantF(DeterminantF);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.D);
    const Vector& r_span = GetValue(KNOT_SPAN_SIZES);
    double squared_length = 0.0;
    for (IndexType i = 0; i < mDim; ++i) squared_length += r_span[i] * r_span[i];
    rValues.SetCharacteristicGeometryLength(0.5 * std::sqrt(squared_length));
}

void SbmSolidGreenLagrangeCondition::CalculateTractionAndTangent(
    const Matrix& rDN_DX, const Matrix& rF, const Matrix& rB,
    const ConstitutiveVariables& rVariables, Vector& rTraction, Matrix& rTractionTangent) const
{
    const SizeType dofs = GetGeometry().size() * mDim;
    const Matrix stress = MathUtils<double>::StressVectorToTensor(rVariables.StressVector);
    const Matrix DB = prod(rVariables.D, rB);
    Vector stress_normal = ZeroVector(mDim);
    for (IndexType j = 0; j < mDim; ++j) {
        for (IndexType k = 0; k < mDim; ++k) stress_normal[j] += stress(j, k) * mNormalPhysicalSpace[k];
    }
    rTraction = prod(rF, stress_normal);
    rTractionTangent = ZeroMatrix(mDim, dofs);
    for (IndexType a = 0; a < GetGeometry().size(); ++a) {
        double geometric = 0.0;
        for (IndexType j = 0; j < mDim; ++j) geometric += rDN_DX(a, j) * stress_normal[j];
        for (IndexType k = 0; k < mDim; ++k) {
            const IndexType col = a * mDim + k;
            const Vector stress_column = column(DB, col);
            const Matrix stress_derivative = MathUtils<double>::StressVectorToTensor(stress_column);
            for (IndexType i = 0; i < mDim; ++i) {
                rTractionTangent(i, col) = i == k ? geometric : 0.0;
                for (IndexType j = 0; j < mDim; ++j) {
                    for (IndexType l = 0; l < mDim; ++l) {
                        rTractionTangent(i, col) += rF(i, j) * stress_derivative(j, l) * mNormalPhysicalSpace[l];
                    }
                }
            }
        }
    }
}

void SbmSolidGreenLagrangeCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    Matrix gradients, F, B;
    CalculateKinematics(gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, MathUtils<double>::Det(F), values, variables, true);
    if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
        mpConstitutiveLaw->InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    Vector traction;
    CalculateTractionAndTangent(gradients, F, B, variables, traction, mInitialTractionTangent);
}

void SbmSolidGreenLagrangeCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType nodes = GetGeometry().size();
    const SizeType dofs = nodes * mDim;
    Matrix gradients, F, B;
    CalculateKinematics(gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, MathUtils<double>::Det(F), values, variables, true);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    Vector traction;
    Matrix tangent;
    CalculateTractionAndTangent(gradients, F, B, variables, traction, tangent);
    KRATOS_ERROR_IF(mInitialTractionTangent.size1() != mDim || mInitialTractionTangent.size2() != dofs)
        << "InitializeSolutionStep must be called before assembling condition " << Id() << std::endl;
    Vector H = ZeroVector(nodes);
    ComputeTaylorExpansionContribution(H);
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
    const double weight = GetValue(INTEGRATION_WEIGHT);
    Matrix projector = IdentityMatrix(mDim);
    if (mpProjectionNode->Has(DIRECTION)) {
        const Vector& r_direction = mpProjectionNode->GetValue(DIRECTION);
        double norm = 0.0;
        for (IndexType i = 0; i < mDim; ++i) norm += r_direction[i] * r_direction[i];
        KRATOS_ERROR_IF(norm <= 1e-24) << "Zero prescribed direction in condition " << Id() << std::endl;
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j) projector(i,j) = r_direction[i] * r_direction[j] / norm;
    }
    rLeftHandSideMatrix = ZeroMatrix(dofs, dofs);
    for (IndexType a = 0; a < nodes; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            for (IndexType b = 0; b < nodes; ++b) {
                for (IndexType j = 0; j < mDim; ++j) {
                    const IndexType col = b * mDim + j;
                    double flux = 0.0, adjoint = 0.0;
                    for (IndexType k = 0; k < mDim; ++k) {
                        flux += projector(i,k) * tangent(k,col);
                        adjoint += mInitialTractionTangent(k,row) * projector(k,j);
                    }
                    rLeftHandSideMatrix(row,col) += weight * (mPenalty * r_N(0,a) * H[b] * projector(i,j)
                        - r_N(0,a) * flux - mNitschePenalty * H[b] * adjoint);
                }
            }
        }
    }
}

void SbmSolidGreenLagrangeCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType nodes = GetGeometry().size();
    const SizeType dofs = nodes * mDim;
    Matrix gradients, F, B;
    CalculateKinematics(gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, MathUtils<double>::Det(F), values, variables, true);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    Vector traction;
    Matrix tangent;
    CalculateTractionAndTangent(gradients, F, B, variables, traction, tangent);
    KRATOS_ERROR_IF(mInitialTractionTangent.size1() != mDim || mInitialTractionTangent.size2() != dofs)
        << "InitializeSolutionStep must be called before assembling condition " << Id() << std::endl;
    Vector H = ZeroVector(nodes);
    ComputeTaylorExpansionContribution(H);
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
    Vector gap = ZeroVector(mDim);
    const auto& r_prescribed = mpProjectionNode->GetValue(DISPLACEMENT);
    for (IndexType i = 0; i < mDim; ++i) gap[i] = r_prescribed[i];
    for (IndexType a = 0; a < nodes; ++a) {
        const auto& r_u = GetGeometry()[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i) gap[i] -= H[a] * r_u[i];
    }
    Matrix projector = IdentityMatrix(mDim);
    if (mpProjectionNode->Has(DIRECTION)) {
        const Vector& r_direction = mpProjectionNode->GetValue(DIRECTION);
        double norm = 0.0;
        for (IndexType i = 0; i < mDim; ++i) norm += r_direction[i] * r_direction[i];
        KRATOS_ERROR_IF(norm <= 1e-24) << "Zero prescribed direction in condition " << Id() << std::endl;
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j) projector(i,j) = r_direction[i] * r_direction[j] / norm;
    }
    const Vector projected_gap = prod(projector, gap);
    const Vector projected_traction = prod(projector, traction);
    const double weight = GetValue(INTEGRATION_WEIGHT);
    rRightHandSideVector = ZeroVector(dofs);
    for (IndexType a = 0; a < nodes; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            double adjoint = 0.0;
            for (IndexType k = 0; k < mDim; ++k) adjoint += mInitialTractionTangent(k,row) * projected_gap[k];
            rRightHandSideVector[row] += weight * (r_N(0,a) * (mPenalty * projected_gap[i] + projected_traction[i])
                - mNitschePenalty * adjoint);
        }
    }
}

void SbmSolidGreenLagrangeCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    Matrix gradients, F, B;
    CalculateKinematics(gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, MathUtils<double>::Det(F), values, variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
    const Vector cauchy_stress = variables.StressVector;
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    const Matrix pk2 = MathUtils<double>::StressVectorToTensor(variables.StressVector);
    Matrix inverse_f(mDim, mDim);
    double det_f = 0.0;
    MathUtils<double>::InvertMatrix(F, inverse_f, det_f);
    array_1d<double, 3> current_normal = ZeroVector(3);
    for (IndexType i = 0; i < mDim; ++i) {
        for (IndexType j = 0; j < mDim; ++j) {
            current_normal[i] += inverse_f(j, i) * mNormalPhysicalSpace[j];
        }
    }
    current_normal /= norm_2(current_normal);
    Vector traction;
    CalculateTraction(cauchy_stress, current_normal, traction);
    SetValue(NORMAL_STRESS, traction);
    SetValue(PK2_STRESS_VECTOR, variables.StressVector);
    SetValue(PK2_STRESS_TENSOR, pk2);
    if (mDim == 2) {
        SetValue(CAUCHY_STRESS_XX, cauchy_stress[0]);
        SetValue(CAUCHY_STRESS_YY, cauchy_stress[1]);
        SetValue(CAUCHY_STRESS_XY, cauchy_stress[2]);
        SetValue(PK2_STRESS_XX, variables.StressVector[0]);
        SetValue(PK2_STRESS_YY, variables.StressVector[1]);
        SetValue(PK2_STRESS_XY, variables.StressVector[2]);
    } else {
        SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(cauchy_stress));
    }
    if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLaw->FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

} // namespace Kratos
