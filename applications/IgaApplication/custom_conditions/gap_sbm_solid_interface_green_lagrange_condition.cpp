// KRATOS Multiphysics - IgaApplication
// License: BSD License

#include <cmath>
#include "custom_conditions/gap_sbm_solid_interface_green_lagrange_condition.h"

namespace Kratos
{

Condition::Pointer GapSbmSolidInterfaceGreenLagrangeCondition::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidInterfaceGreenLagrangeCondition>(NewId, pGeometry, pProperties);
}

Condition::Pointer GapSbmSolidInterfaceGreenLagrangeCondition::Create(
    IndexType NewId, NodesArrayType const& rNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidInterfaceGreenLagrangeCondition>(
        NewId, GetGeometry().Create(rNodes), pProperties);
}

void GapSbmSolidInterfaceGreenLagrangeCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    GapSbmSolidInterfaceCondition::Initialize(rCurrentProcessInfo);
}

void GapSbmSolidInterfaceGreenLagrangeCondition::CalculateSideKinematics(
    const GeometryType& rGeometry, const Vector& rDistance, const Vector& rDisplacement,
    Matrix& rGradients, Matrix& rF, Matrix& rB)
{
    const SizeType nodes = rGeometry.size();
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix transposed_gradients = ZeroMatrix(3, nodes);
    ComputeGradientTaylorExpansionContribution(rGeometry, rDistance, transposed_gradients);
    rGradients = ZeroMatrix(nodes, mDim);
    for (IndexType a = 0; a < nodes; ++a)
        for (IndexType i = 0; i < mDim; ++i)
            rGradients(a, i) = transposed_gradients(i, a);

    rF = IdentityMatrix(mDim);
    for (IndexType a = 0; a < nodes; ++a)
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j)
                rF(i, j) += rDisplacement[a * mDim + i] * rGradients(a, j);
    const double det_f = MathUtils<double>::Det(rF);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_f) && det_f > 0.0)
        << "Invalid deformation gradient determinant " << det_f
        << " in interface condition " << Id() << std::endl;

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

void GapSbmSolidInterfaceGreenLagrangeCondition::SetFiniteStrainParameters(
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
    if (Has(CHARACTERISTIC_GEOMETRY_LENGTH))
        rValues.SetCharacteristicGeometryLength(norm_2(GetValue(CHARACTERISTIC_GEOMETRY_LENGTH)));
}

void GapSbmSolidInterfaceGreenLagrangeCondition::CalculateNominalTractionAndTangent(
    const Matrix& rGradients, const Matrix& rF, const Matrix& rB,
    const ConstitutiveVariables& rVariables, Vector& rTraction, Matrix& rTangent) const
{
    const Matrix stress = MathUtils<double>::StressVectorToTensor(rVariables.StressVector);
    const Matrix DB = prod(rVariables.D, rB);
    Vector normal(mDim);
    for (IndexType i = 0; i < mDim; ++i) normal[i] = mNormalPhysicalSpace[i];
    const Vector stress_normal = prod(stress, normal);
    rTraction = prod(rF, stress_normal);
    rTangent = ZeroMatrix(mDim, rB.size2());
    for (IndexType a = 0; a < rGradients.size1(); ++a) {
        double geometric = 0.0;
        for (IndexType j = 0; j < mDim; ++j) geometric += rGradients(a, j) * stress_normal[j];
        for (IndexType k = 0; k < mDim; ++k) {
            const IndexType col = a * mDim + k;
            const Matrix stress_derivative = MathUtils<double>::StressVectorToTensor(column(DB, col));
            for (IndexType i = 0; i < mDim; ++i) {
                rTangent(i, col) = i == k ? geometric : 0.0;
                for (IndexType j = 0; j < mDim; ++j)
                    for (IndexType l = 0; l < mDim; ++l)
                        rTangent(i, col) += rF(i, j) * stress_derivative(j, l) * normal[l];
            }
        }
    }
}

void GapSbmSolidInterfaceGreenLagrangeCondition::CalculateLeftHandSide(
    MatrixType& rLHS, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_plus_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_minus_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[1];
    const SizeType nodes_plus = r_plus_geometry.size();
    const SizeType nodes_minus = r_minus_geometry.size();
    const SizeType dofs_plus = nodes_plus * mDim;
    const SizeType dofs_minus = nodes_minus * mDim;
    const SizeType shift = dofs_plus;
    Vector H_plus(nodes_plus), H_minus(nodes_minus), u_plus(dofs_plus), u_minus(dofs_minus);
    ComputeTaylorExpansionContribution(r_plus_geometry, mDistanceVectorPlus, H_plus);
    ComputeTaylorExpansionContribution(r_minus_geometry, mDistanceVectorMinus, H_minus);
    GetSolutionCoefficientVectorPlus(u_plus);
    GetSolutionCoefficientVectorMinus(u_minus);
    Matrix gradients_plus, F_plus, B_plus, gradients_minus, F_minus, B_minus;
    CalculateSideKinematics(r_plus_geometry, mDistanceVectorPlus, u_plus, gradients_plus, F_plus, B_plus);
    CalculateSideKinematics(r_minus_geometry, mDistanceVectorMinus, u_minus, gradients_minus, F_minus, B_minus);

    ConstitutiveVariables vars_plus(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveVariables vars_minus(mpConstitutiveLawMinus->GetStrainSize());
    ConstitutiveLaw::Parameters values_plus(r_plus_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveLaw::Parameters values_minus(r_minus_geometry, GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F_plus, values_plus, vars_plus, true);
    SetFiniteStrainParameters(F_minus, values_minus, vars_minus, true);
    mpConstitutiveLaw->CalculateMaterialResponse(values_plus, ConstitutiveLaw::StressMeasure_PK2);
    mpConstitutiveLawMinus->CalculateMaterialResponse(values_minus, ConstitutiveLaw::StressMeasure_PK2);
    Vector traction_plus, traction_minus;
    Matrix tangent_plus, tangent_minus;
    CalculateNominalTractionAndTangent(gradients_plus, F_plus, B_plus, vars_plus, traction_plus, tangent_plus);
    CalculateNominalTractionAndTangent(gradients_minus, F_minus, B_minus, vars_minus, traction_minus, tangent_minus);

    const SizeType total_dofs = dofs_plus + dofs_minus;
    rLHS = ZeroMatrix(total_dofs, total_dofs);
    const double weight = GetValue(INTEGRATION_WEIGHT);
    KRATOS_ERROR_IF(mInitialTractionTangentPlus.size1() != mDim
        || mInitialTractionTangentPlus.size2() != dofs_plus
        || mInitialTractionTangentMinus.size1() != mDim
        || mInitialTractionTangentMinus.size2() != dofs_minus)
        << "InitializeSolutionStep must be called before assembling interface condition " << Id() << std::endl;

    for (IndexType a = 0; a < nodes_plus; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            for (IndexType b = 0; b < nodes_plus; ++b) {
                for (IndexType j = 0; j < mDim; ++j) {
                    const IndexType col = b * mDim + j;
                    rLHS(row, col) += weight * (-0.5 * H_plus[a] * tangent_plus(i, col)
                        + mPenalty * H_plus[a] * H_plus[b] * (i == j)
                        - 0.5 * mNitschePenalty * H_plus[b] * mInitialTractionTangentPlus(j, row));
                }
            }
            for (IndexType b = 0; b < nodes_minus; ++b) {
                for (IndexType j = 0; j < mDim; ++j) {
                    const IndexType col = b * mDim + j;
                    rLHS(row, shift + col) += weight * (-0.5 * H_plus[a] * tangent_minus(i, col)
                        - mPenalty * H_plus[a] * H_minus[b] * (i == j)
                        + 0.5 * mNitschePenalty * H_minus[b] * mInitialTractionTangentPlus(j, row));
                }
            }
        }
    }
    for (IndexType a = 0; a < nodes_minus; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = shift + a * mDim + i;
            const IndexType local_row = a * mDim + i;
            for (IndexType b = 0; b < nodes_plus; ++b) {
                for (IndexType j = 0; j < mDim; ++j) {
                    const IndexType col = b * mDim + j;
                    rLHS(row, col) += weight * (0.5 * H_minus[a] * tangent_plus(i, col)
                        - mPenalty * H_minus[a] * H_plus[b] * (i == j)
                        - 0.5 * mNitschePenalty * H_plus[b] * mInitialTractionTangentMinus(j, local_row));
                }
            }
            for (IndexType b = 0; b < nodes_minus; ++b) {
                for (IndexType j = 0; j < mDim; ++j) {
                    const IndexType col = b * mDim + j;
                    rLHS(row, shift + col) += weight * (0.5 * H_minus[a] * tangent_minus(i, col)
                        + mPenalty * H_minus[a] * H_minus[b] * (i == j)
                        + 0.5 * mNitschePenalty * H_minus[b] * mInitialTractionTangentMinus(j, local_row));
                }
            }
        }
    }
}

void GapSbmSolidInterfaceGreenLagrangeCondition::CalculateRightHandSide(
    VectorType& rRHS, const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_plus_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_minus_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[1];
    const SizeType nodes_plus = r_plus_geometry.size();
    const SizeType nodes_minus = r_minus_geometry.size();
    const SizeType dofs_plus = nodes_plus * mDim;
    const SizeType dofs_minus = nodes_minus * mDim;
    const SizeType shift = dofs_plus;
    Vector H_plus(nodes_plus), H_minus(nodes_minus), u_plus(dofs_plus), u_minus(dofs_minus);
    ComputeTaylorExpansionContribution(r_plus_geometry, mDistanceVectorPlus, H_plus);
    ComputeTaylorExpansionContribution(r_minus_geometry, mDistanceVectorMinus, H_minus);
    GetSolutionCoefficientVectorPlus(u_plus);
    GetSolutionCoefficientVectorMinus(u_minus);
    Matrix gradients_plus, F_plus, B_plus, gradients_minus, F_minus, B_minus;
    CalculateSideKinematics(r_plus_geometry, mDistanceVectorPlus, u_plus, gradients_plus, F_plus, B_plus);
    CalculateSideKinematics(r_minus_geometry, mDistanceVectorMinus, u_minus, gradients_minus, F_minus, B_minus);
    ConstitutiveVariables vars_plus(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveVariables vars_minus(mpConstitutiveLawMinus->GetStrainSize());
    ConstitutiveLaw::Parameters values_plus(r_plus_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveLaw::Parameters values_minus(r_minus_geometry, GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F_plus, values_plus, vars_plus, false);
    SetFiniteStrainParameters(F_minus, values_minus, vars_minus, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values_plus, ConstitutiveLaw::StressMeasure_PK2);
    mpConstitutiveLawMinus->CalculateMaterialResponse(values_minus, ConstitutiveLaw::StressMeasure_PK2);
    Vector traction_plus, traction_minus;
    Matrix tangent_plus, tangent_minus;
    CalculateNominalTractionAndTangent(gradients_plus, F_plus, B_plus, vars_plus, traction_plus, tangent_plus);
    CalculateNominalTractionAndTangent(gradients_minus, F_minus, B_minus, vars_minus, traction_minus, tangent_minus);
    Vector gap = ZeroVector(mDim);
    for (IndexType a = 0; a < nodes_plus; ++a)
        for (IndexType i = 0; i < mDim; ++i) gap[i] += H_plus[a] * u_plus[a * mDim + i];
    for (IndexType a = 0; a < nodes_minus; ++a)
        for (IndexType i = 0; i < mDim; ++i) gap[i] -= H_minus[a] * u_minus[a * mDim + i];

    const Vector average_traction = 0.5 * (traction_plus + traction_minus);
    rRHS = ZeroVector(dofs_plus + dofs_minus);
    const double weight = GetValue(INTEGRATION_WEIGHT);
    for (IndexType a = 0; a < nodes_plus; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = a * mDim + i;
            double adjoint = 0.0;
            for (IndexType j = 0; j < mDim; ++j)
                adjoint += 0.5 * mNitschePenalty * gap[j] * mInitialTractionTangentPlus(j, row);
            rRHS[row] += weight * (H_plus[a] * (average_traction[i] - mPenalty * gap[i]) + adjoint);
        }
    }
    for (IndexType a = 0; a < nodes_minus; ++a) {
        for (IndexType i = 0; i < mDim; ++i) {
            const IndexType row = shift + a * mDim + i;
            const IndexType local_row = a * mDim + i;
            double adjoint = 0.0;
            for (IndexType j = 0; j < mDim; ++j)
                adjoint += 0.5 * mNitschePenalty * gap[j] * mInitialTractionTangentMinus(j, local_row);
            rRHS[row] += weight * (-H_minus[a] * average_traction[i] + mPenalty * H_minus[a] * gap[i] + adjoint);
        }
    }
}

void GapSbmSolidInterfaceGreenLagrangeCondition::UpdateFiniteStrainResponses(
    const ProcessInfo& rCurrentProcessInfo, const bool Finalize)
{
    const auto& r_plus_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_minus_geometry = *GetValue(NEIGHBOUR_GEOMETRIES)[1];
    for (const bool is_minus : {false, true}) {
        const auto& r_geometry = is_minus ? r_minus_geometry : r_plus_geometry;
        const auto& r_distance = is_minus ? mDistanceVectorMinus : mDistanceVectorPlus;
        auto& r_law = is_minus ? *mpConstitutiveLawMinus : *mpConstitutiveLaw;
        Vector displacement(r_geometry.size() * mDim);
        if (is_minus) GetSolutionCoefficientVectorMinus(displacement);
        else GetSolutionCoefficientVectorPlus(displacement);
        Matrix gradients, F, B;
        CalculateSideKinematics(r_geometry, r_distance, displacement, gradients, F, B);
        ConstitutiveVariables variables(r_law.GetStrainSize());
        ConstitutiveLaw::Parameters values(r_geometry, GetProperties(), rCurrentProcessInfo);
        SetFiniteStrainParameters(F, values, variables, !Finalize);
        if (Finalize) {
            r_law.CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
            if (r_law.RequiresFinalizeMaterialResponse()) {
                r_law.FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
            }
        } else {
            r_law.CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
            Vector traction;
            Matrix tangent;
            CalculateNominalTractionAndTangent(gradients, F, B, variables, traction, tangent);
            if (is_minus) mInitialTractionTangentMinus = tangent;
            else mInitialTractionTangentPlus = tangent;
            if (r_law.RequiresInitializeMaterialResponse()) {
                r_law.InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
            }
        }
    }
}

void GapSbmSolidInterfaceGreenLagrangeCondition::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    UpdateFiniteStrainResponses(rCurrentProcessInfo, false);
}

void GapSbmSolidInterfaceGreenLagrangeCondition::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    UpdateFiniteStrainResponses(rCurrentProcessInfo, true);
}

} // namespace Kratos
