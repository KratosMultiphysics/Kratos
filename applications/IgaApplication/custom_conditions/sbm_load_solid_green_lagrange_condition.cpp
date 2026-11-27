// KRATOS Multiphysics - IgaApplication
// License: BSD License

#include <cmath>
#include "custom_conditions/sbm_load_solid_green_lagrange_condition.h"

namespace Kratos
{

Condition::Pointer SbmLoadSolidGreenLagrangeCondition::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmLoadSolidGreenLagrangeCondition>(NewId, pGeometry, pProperties);
}

Condition::Pointer SbmLoadSolidGreenLagrangeCondition::Create(
    IndexType NewId, NodesArrayType const& rNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmLoadSolidGreenLagrangeCondition>(
        NewId, GetGeometry().Create(rNodes), pProperties);
}

void SbmLoadSolidGreenLagrangeCondition::CalculateReferenceGradients(
    Matrix& rSurrogateGradients, Matrix& rProjectedGradients)
{
    const auto& r_geometry = GetGeometry();
    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "Unsupported dimension " << mDim << std::endl;
    const auto& r_local_gradients = r_geometry.ShapeFunctionsLocalGradients(GetIntegrationMethod());
    Matrix jacobian;
    jacobian = IdentityMatrix(mDim); // DEBUG: assume an identity parametric-to-physical mapping.
    Matrix inverse_jacobian(mDim, mDim);
    double determinant = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, inverse_jacobian, determinant);
    rSurrogateGradients = prod(r_local_gradients[0], inverse_jacobian);

    Matrix taylor_gradient(3, r_geometry.size());
    ComputeGradientTaylorExpansionContribution(mDistanceVector, taylor_gradient);
    Matrix projected_local(r_geometry.size(), mDim);
    for (IndexType a = 0; a < r_geometry.size(); ++a) {
        for (IndexType i = 0; i < mDim; ++i) projected_local(a, i) = taylor_gradient(i, a);
    }
    rProjectedGradients = prod(projected_local, inverse_jacobian);
}

void SbmLoadSolidGreenLagrangeCondition::CalculateDeformationGradientAndB(
    const Matrix& rGradients, Matrix& rF, Matrix& rB) const
{
    const SizeType nodes = GetGeometry().size();
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    KRATOS_ERROR_IF(strain_size != (mDim == 2 ? 3 : 6)) << "Incompatible constitutive strain size" << std::endl;
    rF = IdentityMatrix(mDim);
    for (IndexType a = 0; a < nodes; ++a) {
        const auto& r_displacement = GetGeometry()[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) rF(i, j) += r_displacement[i] * rGradients(a, j);
        }
    }
    const double determinant = MathUtils<double>::Det(rF);
    KRATOS_ERROR_IF_NOT(std::isfinite(determinant) && determinant > 0.0)
        << "Invalid deformation gradient determinant " << determinant << " in condition " << Id() << std::endl;
    rB = ZeroMatrix(strain_size, nodes * mDim);
    for (IndexType a = 0; a < nodes; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType col = a * mDim + i;
            rB(0, col) = rF(i, 0) * rGradients(a, 0);
            rB(1, col) = rF(i, 1) * rGradients(a, 1);
            if (mDim == 2) {
                rB(2, col) = rF(i, 0) * rGradients(a, 1) + rF(i, 1) * rGradients(a, 0);
            } else {
                rB(2, col) = rF(i, 2) * rGradients(a, 2);
                rB(3, col) = rF(i, 0) * rGradients(a, 1) + rF(i, 1) * rGradients(a, 0);
                rB(4, col) = rF(i, 1) * rGradients(a, 2) + rF(i, 2) * rGradients(a, 1);
                rB(5, col) = rF(i, 0) * rGradients(a, 2) + rF(i, 2) * rGradients(a, 0);
            }
        }
    }
}

void SbmLoadSolidGreenLagrangeCondition::SetConstitutiveParameters(
    const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
    ConstitutiveVariables& rVariables, bool ComputeTangent) const
{
    auto& r_options = rValues.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, ComputeTangent);
    rValues.SetDeformationGradientF(rF);
    rValues.SetDeterminantF(MathUtils<double>::Det(rF));
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.D);
    const Vector& r_span = GetValue(KNOT_SPAN_SIZES);
    double squared_length = 0.0;
    for (IndexType i = 0; i < mDim; ++i) squared_length += r_span[i] * r_span[i];
    rValues.SetCharacteristicGeometryLength(0.5 * std::sqrt(squared_length));
}

void SbmLoadSolidGreenLagrangeCondition::CalculateNominalTraction(
    const Matrix& rF, const Vector& rPK2Stress,
    const array_1d<double, 3>& rNormal, Vector& rTraction) const
{
    const Matrix stress = MathUtils<double>::StressVectorToTensor(rPK2Stress);
    Vector stress_normal = ZeroVector(mDim);
    for (IndexType i = 0; i < mDim; ++i) {
        for (IndexType j = 0; j < mDim; ++j) stress_normal[i] += stress(i, j) * rNormal[j];
    }
    rTraction = prod(rF, stress_normal);
}

void SbmLoadSolidGreenLagrangeCondition::CalculateNominalTractionTangent(
    const Matrix& rGradients, const Matrix& rF, const Matrix& rB,
    const ConstitutiveVariables& rVariables, const array_1d<double, 3>& rNormal,
    Matrix& rTangent) const
{
    const SizeType dofs = GetGeometry().size() * mDim;
    const Matrix stress = MathUtils<double>::StressVectorToTensor(rVariables.StressVector);
    const Matrix DB = prod(rVariables.D, rB);
    Vector stress_normal = ZeroVector(mDim);
    for (IndexType j = 0; j < mDim; ++j) {
        for (IndexType k = 0; k < mDim; ++k) stress_normal[j] += stress(j, k) * rNormal[k];
    }
    rTangent = ZeroMatrix(mDim, dofs);
    for (IndexType a = 0; a < GetGeometry().size(); ++a) {
        double geometric = 0.0;
        for (IndexType j = 0; j < mDim; ++j) geometric += rGradients(a, j) * stress_normal[j];
        for (IndexType k = 0; k < mDim; ++k) {
            const IndexType col = a * mDim + k;
            const Vector stress_column = column(DB, col);
            const Matrix stress_derivative = MathUtils<double>::StressVectorToTensor(stress_column);
            for (IndexType i = 0; i < mDim; ++i) {
                rTangent(i, col) = i == k ? geometric : 0.0;
                for (IndexType j = 0; j < mDim; ++j) {
                    for (IndexType l = 0; l < mDim; ++l)
                        rTangent(i, col) += rF(i, j) * stress_derivative(j, l) * rNormal[l];
                }
            }
        }
    }
}

void SbmLoadSolidGreenLagrangeCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dofs = r_geometry.size() * mDim;
    Matrix gradients, projected_gradients, F, projected_F, B, projected_B;
    CalculateReferenceGradients(gradients, projected_gradients);
    CalculateDeformationGradientAndB(gradients, F, B);
    CalculateDeformationGradientAndB(projected_gradients, projected_F, projected_B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveVariables projected_variables(mpConstitutiveLawAtProjectionPoint->GetStrainSize());
    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveLaw::Parameters projected_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, true);
    SetConstitutiveParameters(projected_F, projected_values, projected_variables, true);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    mpConstitutiveLawAtProjectionPoint->CalculateMaterialResponse(projected_values, ConstitutiveLaw::StressMeasure_PK2);
    Matrix tangent, projected_tangent;
    CalculateNominalTractionTangent(gradients, F, B, variables, mNormalPhysicalSpace, tangent);
    CalculateNominalTractionTangent(projected_gradients, projected_F, projected_B, projected_variables, mTrueNormal, projected_tangent);
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());
    const double weight = GetValue(INTEGRATION_WEIGHT);
    rLeftHandSideMatrix = ZeroMatrix(dofs, dofs);
    for (IndexType a = 0; a < r_geometry.size(); ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            for (IndexType col = 0; col < dofs; ++col) {
                rLeftHandSideMatrix(row, col) += weight * r_N(0, a) *
                    (-tangent(i, col) + mTrueDotSurrogateNormal * projected_tangent(i, col));
            }
        }
    }
}

void SbmLoadSolidGreenLagrangeCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dofs = r_geometry.size() * mDim;
    Matrix gradients, projected_gradients, F, projected_F, B, projected_B;
    CalculateReferenceGradients(gradients, projected_gradients);
    CalculateDeformationGradientAndB(gradients, F, B);
    CalculateDeformationGradientAndB(projected_gradients, projected_F, projected_B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveVariables projected_variables(mpConstitutiveLawAtProjectionPoint->GetStrainSize());
    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveLaw::Parameters projected_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, false);
    SetConstitutiveParameters(projected_F, projected_values, projected_variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    mpConstitutiveLawAtProjectionPoint->CalculateMaterialResponse(projected_values, ConstitutiveLaw::StressMeasure_PK2);
    Vector traction, projected_traction;
    CalculateNominalTraction(F, variables.StressVector, mNormalPhysicalSpace, traction);
    CalculateNominalTraction(projected_F, projected_variables.StressVector, mTrueNormal, projected_traction);
    const auto& r_force = mpProjectionNode->GetValue(FORCE);
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());
    const double weight = GetValue(INTEGRATION_WEIGHT);
    rRightHandSideVector = ZeroVector(dofs);
    for (IndexType a = 0; a < r_geometry.size(); ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            rRightHandSideVector[a * mDim + i] += weight * r_N(0, a) *
                (traction[i] + mTrueDotSurrogateNormal * (r_force[i] - projected_traction[i]));
        }
    }
}

void SbmLoadSolidGreenLagrangeCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    Matrix gradients, projected_gradients, F, projected_F, B, projected_B;
    CalculateReferenceGradients(gradients, projected_gradients);
    CalculateDeformationGradientAndB(gradients, F, B);
    CalculateDeformationGradientAndB(projected_gradients, projected_F, projected_B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveVariables projected_variables(mpConstitutiveLawAtProjectionPoint->GetStrainSize());
    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveLaw::Parameters projected_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, false);
    SetConstitutiveParameters(projected_F, projected_values, projected_variables, false);
    if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
        mpConstitutiveLaw->InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
    if (mpConstitutiveLawAtProjectionPoint->RequiresInitializeMaterialResponse()) {
        mpConstitutiveLawAtProjectionPoint->InitializeMaterialResponse(projected_values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

void SbmLoadSolidGreenLagrangeCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    Matrix gradients, projected_gradients, F, projected_F, B, projected_B;
    CalculateReferenceGradients(gradients, projected_gradients);
    CalculateDeformationGradientAndB(gradients, F, B);
    CalculateDeformationGradientAndB(projected_gradients, projected_F, projected_B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveVariables projected_variables(mpConstitutiveLawAtProjectionPoint->GetStrainSize());
    ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveLaw::Parameters projected_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, false);
    SetConstitutiveParameters(projected_F, projected_values, projected_variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
    const Vector cauchy = variables.StressVector;
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    mpConstitutiveLawAtProjectionPoint->CalculateMaterialResponse(projected_values, ConstitutiveLaw::StressMeasure_PK2);
    SetValue(PK2_STRESS_VECTOR, variables.StressVector);
    SetValue(PK2_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(variables.StressVector));
    Matrix inverse_f(mDim, mDim);
    double determinant = 0.0;
    MathUtils<double>::InvertMatrix(F, inverse_f, determinant);
    array_1d<double, 3> current_normal = ZeroVector(3);
    for (IndexType i = 0; i < mDim; ++i)
        for (IndexType j = 0; j < mDim; ++j) current_normal[i] += inverse_f(j, i) * mNormalPhysicalSpace[j];
    current_normal /= norm_2(current_normal);
    Vector traction;
    CalculateTraction(cauchy, current_normal, traction);
    SetValue(NORMAL_STRESS, traction);
    if (mDim == 2) {
        SetValue(CAUCHY_STRESS_XX, cauchy[0]);
        SetValue(CAUCHY_STRESS_YY, cauchy[1]);
        SetValue(CAUCHY_STRESS_XY, cauchy[2]);
        SetValue(PK2_STRESS_XX, variables.StressVector[0]);
        SetValue(PK2_STRESS_YY, variables.StressVector[1]);
        SetValue(PK2_STRESS_XY, variables.StressVector[2]);
    } else {
        SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(cauchy));
    }
    if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLaw->FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
    if (mpConstitutiveLawAtProjectionPoint->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLawAtProjectionPoint->FinalizeMaterialResponse(projected_values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

} // namespace Kratos
