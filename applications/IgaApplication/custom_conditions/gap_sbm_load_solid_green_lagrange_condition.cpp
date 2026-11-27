// KRATOS Multiphysics - IgaApplication
// License: BSD License

#include <cmath>
#include "custom_conditions/gap_sbm_load_solid_green_lagrange_condition.h"

namespace Kratos
{

Condition::Pointer GapSbmLoadSolidGreenLagrangeCondition::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmLoadSolidGreenLagrangeCondition>(NewId, pGeometry, pProperties);
}

Condition::Pointer GapSbmLoadSolidGreenLagrangeCondition::Create(
    IndexType NewId, NodesArrayType const& rNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmLoadSolidGreenLagrangeCondition>(
        NewId, GetGeometry().Create(rNodes), pProperties);
}

void GapSbmLoadSolidGreenLagrangeCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_GEOMETRIES) && !GetValue(NEIGHBOUR_GEOMETRIES).empty())
        << "GapSbmLoadSolidGreenLagrangeCondition " << Id()
        << " requires a surrogate NEIGHBOUR_GEOMETRIES entry; check which modeler created the condition."
        << std::endl;
    GapSbmLoadSolidCondition::Initialize(rCurrentProcessInfo);
}

void GapSbmLoadSolidGreenLagrangeCondition::CalculateGreenKinematics(
    Matrix& rGradients, Matrix& rF, Matrix& rB)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType nodes = r_geometry.size();
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix transposed_gradients = ZeroMatrix(3, nodes);
    ComputeGradientTaylorExpansionContribution(transposed_gradients);
    rGradients = ZeroMatrix(nodes, mDim);
    for (IndexType a = 0; a < nodes; ++a)
        for (IndexType i = 0; i < mDim; ++i)
            rGradients(a, i) = transposed_gradients(i, a);

    rF = IdentityMatrix(mDim);
    for (IndexType a = 0; a < nodes; ++a) {
        const auto& r_displacement = r_geometry[a].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType i = 0; i < mDim; ++i)
            for (IndexType j = 0; j < mDim; ++j)
                rF(i, j) += r_displacement[i] * rGradients(a, j);
    }
    const double determinant_f = MathUtils<double>::Det(rF);
    KRATOS_ERROR_IF_NOT(std::isfinite(determinant_f) && determinant_f > 0.0)
        << "Invalid deformation gradient determinant " << determinant_f
        << " in condition " << Id() << std::endl;

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

void GapSbmLoadSolidGreenLagrangeCondition::SetFiniteStrainParameters(
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

void GapSbmLoadSolidGreenLagrangeCondition::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    Matrix gradients, F, B;
    CalculateGreenKinematics(gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F, values, variables, true);
    if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
        mpConstitutiveLaw->InitializeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

void GapSbmLoadSolidGreenLagrangeCondition::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    Matrix gradients, F, B;
    CalculateGreenKinematics(gradients, F, B);
    ConstitutiveVariables variables(mpConstitutiveLaw->GetStrainSize());
    ConstitutiveLaw::Parameters values(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);
    SetFiniteStrainParameters(F, values, variables, false);
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_Cauchy);
    const Vector cauchy_stress = variables.StressVector;
    mpConstitutiveLaw->CalculateMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    const Vector pk2_stress = variables.StressVector;

    SetValue(PK2_STRESS_VECTOR, pk2_stress);
    SetValue(PK2_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(pk2_stress));
    if (mDim == 2) {
        SetValue(CAUCHY_STRESS_XX, cauchy_stress[0]);
        SetValue(CAUCHY_STRESS_YY, cauchy_stress[1]);
        SetValue(CAUCHY_STRESS_XY, cauchy_stress[2]);
    } else {
        SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(cauchy_stress));
    }

    Matrix inverse_f(mDim, mDim);
    double determinant_f = 0.0;
    MathUtils<double>::InvertMatrix(F, inverse_f, determinant_f);
    Vector reference_normal = ZeroVector(mDim);
    for (IndexType i = 0; i < mDim; ++i) reference_normal[i] = mNormalPhysicalSpace[i];
    Vector current_normal = prod(trans(inverse_f), reference_normal);
    current_normal /= norm_2(current_normal);
    const Matrix cauchy_tensor = MathUtils<double>::StressVectorToTensor(cauchy_stress);
    SetValue(NORMAL_STRESS, prod(cauchy_tensor, current_normal));

    Vector H = ZeroVector(GetSurrogateGeometry().size());
    ComputeTaylorExpansionContribution(H);
    Vector displacement = ZeroVector(3);
    Vector coefficients(GetSurrogateGeometry().size() * mDim);
    GetSolutionCoefficientVector(coefficients);
    for (IndexType a = 0; a < H.size(); ++a)
        for (IndexType i = 0; i < mDim; ++i)
            displacement[i] += H[a] * coefficients[a * mDim + i];
    SetValue(DISPLACEMENT, displacement);

    if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
        mpConstitutiveLaw->FinalizeMaterialResponse(values, ConstitutiveLaw::StressMeasure_PK2);
    }
}

} // namespace Kratos
