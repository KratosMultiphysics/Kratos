// KRATOS Multiphysics - IgaApplication
// License: BSD License

#include <cmath>
#include "custom_conditions/gap_sbm_solid_green_lagrange_condition.h"

namespace Kratos
{

Condition::Pointer GapSbmSolidGreenLagrangeCondition::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidGreenLagrangeCondition>(NewId, pGeometry, pProperties);
}

Condition::Pointer GapSbmSolidGreenLagrangeCondition::Create(
    IndexType NewId, NodesArrayType const& rNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidGreenLagrangeCondition>(
        NewId, GetGeometry().Create(rNodes), pProperties);
}

void GapSbmSolidGreenLagrangeCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_GEOMETRIES) && !GetValue(NEIGHBOUR_GEOMETRIES).empty())
        << "GapSbmSolidGreenLagrangeCondition " << Id()
        << " requires a surrogate NEIGHBOUR_GEOMETRIES entry; check which modeler created the condition."
        << std::endl;
    GapSbmSolidCondition::Initialize(rCurrentProcessInfo);
}

void GapSbmSolidGreenLagrangeCondition::CalculateGapKinematics(
    Vector& rTaylorShapeFunctions, Matrix& rGradients, Matrix& rF, Matrix& rB)
{
    const auto& r_surrogate_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const SizeType nodes = r_surrogate_geometry.size();
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();

    rTaylorShapeFunctions = ZeroVector(nodes);
    ComputeTaylorExpansionContribution(rTaylorShapeFunctions);
    Matrix transposed_gradients = ZeroMatrix(3, nodes);
    ComputeGradientTaylorExpansionContribution(transposed_gradients);
    rGradients = ZeroMatrix(nodes, mDim);
    for (IndexType a = 0; a < nodes; ++a)
        for (IndexType i = 0; i < mDim; ++i)
            rGradients(a, i) = transposed_gradients(i, a);

    rF = IdentityMatrix(mDim);
    for (IndexType a = 0; a < nodes; ++a) {
        const auto& r_displacement = r_surrogate_geometry[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j)
                rF(i, j) += r_displacement[i] * rGradients(a, j);
    }
    const double determinant_f = MathUtils<double>::Det(rF);
    KRATOS_ERROR_IF_NOT(std::isfinite(determinant_f) && determinant_f > 0.0)
        << "Invalid deformation gradient determinant " << determinant_f
        << " in condition " << Id() << std::endl;

    const SizeType dofs = nodes * mDim;
    rB = ZeroMatrix(strain_size, dofs);
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

void GapSbmSolidGreenLagrangeCondition::SetConstitutiveParameters(
    const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
    ConstitutiveVariables& rVariables, const bool ComputeTangent) const
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
    rValues.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());
}

void GapSbmSolidGreenLagrangeCondition::CalculateNominalTractionTangent(
    const Matrix& rGradients, const Matrix& rF, const Matrix& rB,
    const ConstitutiveVariables& rVariables, Matrix& rTangent) const
{
    const SizeType dofs = rB.size2();
    const Matrix stress = MathUtils<double>::StressVectorToTensor(rVariables.StressVector);
    const Matrix DB = prod(rVariables.D, rB);
    Vector stress_normal = ZeroVector(mDim);
    for (IndexType i = 0; i < mDim; ++i)
        for (IndexType j = 0; j < mDim; ++j)
            stress_normal[i] += stress(i, j) * mNormalPhysicalSpace[j];

    rTangent = ZeroMatrix(mDim, dofs);
    for (IndexType a = 0; a < rGradients.size1(); ++a) {
        double geometric = 0.0;
        for (IndexType j = 0; j < mDim; ++j)
            geometric += rGradients(a, j) * stress_normal[j];
        for (IndexType k = 0; k < mDim; ++k) {
            const IndexType col = a * mDim + k;
            const Matrix stress_derivative = MathUtils<double>::StressVectorToTensor(column(DB, col));
            for (IndexType i = 0; i < mDim; ++i) {
                rTangent(i, col) = (i == k) ? geometric : 0.0;
                for (IndexType j = 0; j < mDim; ++j)
                    for (IndexType l = 0; l < mDim; ++l)
                        rTangent(i, col) += rF(i, j) * stress_derivative(j, l) * mNormalPhysicalSpace[l];
            }
        }
    }
}

void GapSbmSolidGreenLagrangeCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGapKinematics(H, gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, true);
    if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
        mpConstitutiveLaw->InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    CalculateNominalTractionTangent(gradients, F, B, variables, mInitialTractionTangent);
}

void GapSbmSolidGreenLagrangeCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGapKinematics(H, gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, true);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    Matrix tangent;
    CalculateNominalTractionTangent(gradients, F, B, variables, tangent);

    const SizeType nodes = H.size();
    const SizeType dofs = nodes * mDim;
    rLeftHandSideMatrix = ZeroMatrix(dofs, dofs);
    KRATOS_ERROR_IF(mInitialTractionTangent.size1() != mDim || mInitialTractionTangent.size2() != dofs)
        << "InitializeSolutionStep must be called before assembling condition " << Id() << std::endl;

    Matrix projector = IdentityMatrix(mDim);
    if (Has(DIRECTION)) {
        const Vector& r_direction = GetValue(DIRECTION);
        double norm_squared = 0.0;
        for (IndexType i = 0; i < mDim; ++i) norm_squared += r_direction[i] * r_direction[i];
        KRATOS_ERROR_IF(norm_squared <= 1e-24) << "Zero prescribed direction in condition " << Id() << std::endl;
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j)
                projector(i, j) = r_direction[i] * r_direction[j] / norm_squared;
    }

    const double weight = GetValue(INTEGRATION_WEIGHT);
    for (IndexType a = 0; a < nodes; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            for (IndexType b = 0; b < nodes; ++b) {
                for (IndexType j = 0; j < mDim; ++j) {
                    const IndexType col = b * mDim + j;
                    double flux = 0.0, adjoint = 0.0;
                    for (IndexType k = 0; k < mDim; ++k) {
                        flux += projector(i, k) * tangent(k, col);
                        adjoint += mInitialTractionTangent(k, row) * projector(k, j);
                    }
                    rLeftHandSideMatrix(row, col) += weight * (mPenalty * H[a] * H[b] * projector(i, j)
                        - H[a] * flux - mNitschePenalty * H[b] * adjoint);
                }
            }
        }
    }
}

void GapSbmSolidGreenLagrangeCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGapKinematics(H, gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    const Matrix stress = MathUtils<double>::StressVectorToTensor(variables.StressVector);
    Vector normal = ZeroVector(mDim);
    for (IndexType i = 0; i < mDim; ++i) normal[i] = mNormalPhysicalSpace[i];
    const Vector stress_normal = prod(stress, normal);
    const Vector traction = prod(F, stress_normal);
    const SizeType nodes = H.size();
    const SizeType dofs = nodes * mDim;
    KRATOS_ERROR_IF(mInitialTractionTangent.size1() != mDim || mInitialTractionTangent.size2() != dofs)
        << "InitializeSolutionStep must be called before assembling condition " << Id() << std::endl;

    Vector gap = ZeroVector(mDim);
    const Vector& r_prescribed = GetValue(DISPLACEMENT);
    for (IndexType i = 0; i < mDim; ++i) gap[i] = r_prescribed[i];
    const auto& r_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[0];
    for (IndexType a = 0; a < nodes; ++a) {
        const auto& r_u = r_geometry[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i) gap[i] -= H[a] * r_u[i];
    }

    Matrix projector = IdentityMatrix(mDim);
    if (Has(DIRECTION)) {
        const Vector& r_direction = GetValue(DIRECTION);
        double norm_squared = 0.0;
        for (IndexType i = 0; i < mDim; ++i) norm_squared += r_direction[i] * r_direction[i];
        KRATOS_ERROR_IF(norm_squared <= 1e-24) << "Zero prescribed direction in condition " << Id() << std::endl;
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j)
                projector(i, j) = r_direction[i] * r_direction[j] / norm_squared;
    }

    const Vector projected_gap = prod(projector, gap);
    const Vector projected_traction = prod(projector, traction);
    rRightHandSideVector = ZeroVector(dofs);
    const double weight = GetValue(INTEGRATION_WEIGHT);
    for (IndexType a = 0; a < nodes; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            double adjoint = 0.0;
            for (IndexType k = 0; k < mDim; ++k)
                adjoint += mInitialTractionTangent(k, row) * projected_gap[k];
            rRightHandSideVector[row] = weight * (H[a] * (mPenalty * projected_gap[i] + projected_traction[i])
                - mNitschePenalty * adjoint);
        }
    }
}

void GapSbmSolidGreenLagrangeCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGapKinematics(H, gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    SetConstitutiveParameters(F, values, variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
    const Vector cauchy_stress = variables.StressVector;
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    SetValue(PK2_STRESS_VECTOR, variables.StressVector);
    SetValue(PK2_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(variables.StressVector));
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
    Matrix inverse_f(mDim, mDim);
    double determinant_f = 0.0;
    MathUtils<double>::InvertMatrix(F, inverse_f, determinant_f);
    array_1d<double, 3> current_normal = ZeroVector(3);
    for (IndexType i = 0; i < mDim; ++i)
        for (IndexType j = 0; j < mDim; ++j)
            current_normal[i] += inverse_f(j, i) * mNormalPhysicalSpace[j];
    current_normal /= norm_2(current_normal);
    Vector normal_traction = ZeroVector(mDim);
    const Matrix cauchy_tensor = MathUtils<double>::StressVectorToTensor(cauchy_stress);
    for (IndexType i = 0; i < mDim; ++i)
        for (IndexType j = 0; j < mDim; ++j)
            normal_traction[i] += cauchy_tensor(i, j) * current_normal[j];
    SetValue(NORMAL_STRESS, normal_traction);
    if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLaw->FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

} // namespace Kratos
