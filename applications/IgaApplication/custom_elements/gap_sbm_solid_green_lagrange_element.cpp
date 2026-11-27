// KRATOS Multiphysics - IgaApplication
// License: BSD License

#include <cmath>
#include "custom_elements/gap_sbm_solid_green_lagrange_element.h"

namespace Kratos
{

Element::Pointer GapSbmSolidGreenLagrangeElement::Create(
    IndexType NewId, NodesArrayType const& rNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidGreenLagrangeElement>(
        NewId, GetSurrogateGeometry().Create(rNodes), pProperties);
}

Element::Pointer GapSbmSolidGreenLagrangeElement::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidGreenLagrangeElement>(NewId, pGeometry, pProperties);
}

void GapSbmSolidGreenLagrangeElement::CalculateGreenKinematics(
    Vector& rTaylorShapeFunctions, Matrix& rGradients, Matrix& rF, Matrix& rB)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType nodes = r_geometry.size();
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
        const auto& r_u = r_geometry[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j)
                rF(i, j) += r_u[i] * rGradients(a, j);
    }
    const double det_f = MathUtils<double>::Det(rF);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_f) && det_f > 0.0)
        << "Invalid deformation gradient determinant " << det_f
        << " in element " << Id() << std::endl;

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

void GapSbmSolidGreenLagrangeElement::SetFiniteStrainParameters(
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

void GapSbmSolidGreenLagrangeElement::CalculateAndAddGeometricStiffness(
    MatrixType& rLHS, const Matrix& rGradients, const Vector& rPK2Stress,
    const double IntegrationWeight) const
{
    const Matrix weighted_stress = IntegrationWeight * MathUtils<double>::StressVectorToTensor(rPK2Stress);
    Matrix reduced_stiffness(rGradients.size1(), rGradients.size1());
    MathUtils<double>::BDBtProductOperation(reduced_stiffness, weighted_stress, rGradients);
    MathUtils<double>::ExpandAndAddReducedMatrix(rLHS, reduced_stiffness, mDim);
}

void GapSbmSolidGreenLagrangeElement::CalculateLeftHandSide(
    MatrixType& rLHS, const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGreenKinematics(H, gradients, F, B);
    const SizeType dofs = B.size2();
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F, values, variables, true);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);

    rLHS = ZeroMatrix(dofs, dofs);
    const double weight = GetValue(INTEGRATION_WEIGHT);
    noalias(rLHS) += weight * prod(trans(B), Matrix(prod(values.GetConstitutiveMatrix(), B)));
    CalculateAndAddGeometricStiffness(rLHS, gradients, values.GetStressVector(), weight);
}

void GapSbmSolidGreenLagrangeElement::CalculateRightHandSide(
    VectorType& rRHS, const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGreenKinematics(H, gradients, F, B);
    const SizeType nodes = H.size();
    const SizeType dofs = B.size2();
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F, values, variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);

    rRHS = ZeroVector(dofs);
    const double weight = GetValue(INTEGRATION_WEIGHT);
    const Vector body_force = GetValue(BODY_FORCE);
    for (IndexType a = 0; a < nodes; ++a)
        for (IndexType i = 0; i < mDim; ++i)
            rRHS[a * mDim + i] += weight * H[a] * body_force[i];
    noalias(rRHS) -= weight * prod(trans(B), values.GetStressVector());
}

void GapSbmSolidGreenLagrangeElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGreenKinematics(H, gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F, values, variables, true);
    if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
        mpConstitutiveLaw->InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

void GapSbmSolidGreenLagrangeElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    Vector H;
    Matrix gradients, F, B;
    CalculateGreenKinematics(H, gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F, values, variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
    const Vector cauchy = variables.StressVector;
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    SetValue(PK2_STRESS_VECTOR, variables.StressVector);
    SetValue(PK2_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(variables.StressVector));
    if (mDim == 2) {
        SetValue(CAUCHY_STRESS_XX, cauchy[0]);
        SetValue(CAUCHY_STRESS_YY, cauchy[1]);
        SetValue(CAUCHY_STRESS_XY, cauchy[2]);
    } else {
        SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(cauchy));
    }
    if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLaw->FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }

    array_1d<double, 3> displacement = ZeroVector(3);
    Vector coefficients(GetSurrogateGeometry().size() * mDim);
    GetSolutionCoefficientVector(coefficients);
    for (IndexType a = 0; a < H.size(); ++a)
        for (IndexType i = 0; i < mDim; ++i)
            displacement[i] += H[a] * coefficients[a * mDim + i];
    SetValue(DISPLACEMENT, displacement);
}

} // namespace Kratos
