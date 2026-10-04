// License: BSD License
// Kratos default license: kratos/license.txt

#include "small_strain_swift_j2_plasticity_3d.h"

#include <cmath>

#include "constitutive_laws_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainSwiftJ2Plasticity3D::Clone() const
{
    return Kratos::make_shared<SmallStrainSwiftJ2Plasticity3D>(*this);
}

void SmallStrainSwiftJ2Plasticity3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures = {StrainMeasure_Infinitesimal};
    rFeatures.mStrainSize = 6;
    rFeatures.mSpaceDimension = 3;
}

void SmallStrainSwiftJ2Plasticity3D::InitializeMaterial(
    const Properties& rProperties, const GeometryType& rGeometry,
    const Vector& rShapeFunctionsValues)
{
    CheckMaterialProperties(rProperties);
    mPlasticStrain = ZeroVector(6);
    mAccumulatedPlasticStrain = 0.0;
}

void SmallStrainSwiftJ2Plasticity3D::CheckMaterialProperties(const Properties& rProperties)
{
    for (const auto* p_variable : {&YOUNG_MODULUS, &SWIFT_COEFFICIENT,
            &SWIFT_INITIAL_STRAIN, &SWIFT_HARDENING_EXPONENT}) {
        KRATOS_ERROR_IF_NOT(rProperties.Has(*p_variable)) << "Swift J2: missing " << p_variable->Name();
        const double value = rProperties[*p_variable];
        KRATOS_ERROR_IF(!std::isfinite(value) || value <= 0.0)
            << "Swift J2: " << p_variable->Name() << " must be finite and positive; got " << value;
    }
    KRATOS_ERROR_IF_NOT(rProperties.Has(POISSON_RATIO)) << "Swift J2: missing POISSON_RATIO";
    const double nu = rProperties[POISSON_RATIO];
    KRATOS_ERROR_IF(!std::isfinite(nu) || nu <= -1.0 || nu >= 0.5)
        << "Swift J2: POISSON_RATIO must satisfy -1 < nu < 0.5; got " << nu;
    const double young = rProperties[YOUNG_MODULUS];
    KRATOS_ERROR_IF(!std::isfinite(young / (1.0 + nu)) ||
        !std::isfinite(young / (1.0 - 2.0 * nu)) ||
        !std::isfinite(young / ((1.0 + nu) * (1.0 - 2.0 * nu))))
        << "Swift J2: elastic moduli overflow";
    double yield_stress = 0.0;
    double hardening = 0.0;
    EvaluateHardening(rProperties, 0.0, yield_stress, hardening);
}

int SmallStrainSwiftJ2Plasticity3D::Check(const Properties& rProperties,
    const GeometryType& rGeometry, const ProcessInfo& rProcessInfo) const
{
    CheckMaterialProperties(rProperties);
    KRATOS_ERROR_IF(rGeometry.WorkingSpaceDimension() != 3) << "Swift J2 requires a 3D geometry";
    return 0;
}

void SmallStrainSwiftJ2Plasticity3D::EvaluateHardening(const Properties& rProperties,
    const double P, double& rYieldStress, double& rHardeningModulus)
{
    const double base = rProperties[SWIFT_INITIAL_STRAIN] + P;
    KRATOS_ERROR_IF(!std::isfinite(P) || P < 0.0 || !std::isfinite(base) || base <= 0.0)
        << "Swift J2: invalid accumulated plastic strain or power base; p = " << P;
    const double exponent = rProperties[SWIFT_HARDENING_EXPONENT];
    rYieldStress = rProperties[SWIFT_COEFFICIENT] * std::pow(base, exponent);
    rHardeningModulus = exponent * (rYieldStress / base);
    KRATOS_ERROR_IF(!std::isfinite(rYieldStress) || rYieldStress <= 0.0 ||
        !std::isfinite(rHardeningModulus) || rHardeningModulus <= 0.0)
        << "Swift J2: flow stress or hardening modulus is not representable; p = " << P;
}

void SmallStrainSwiftJ2Plasticity3D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    Vector plastic_strain = ZeroVector(6);
    double p = 0.0;
    IntegrateResponse(rValues, plastic_strain, p);
}

void SmallStrainSwiftJ2Plasticity3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    Vector plastic_strain = ZeroVector(6);
    double p = 0.0;
    IntegrateResponse(rValues, plastic_strain, p);
    mPlasticStrain = plastic_strain;
    mAccumulatedPlasticStrain = p;
}

void SmallStrainSwiftJ2Plasticity3D::IntegrateResponse(
    Parameters& rValues, Vector& rPlasticStrain, double& rP)
{
    KRATOS_TRY
    const auto& r_properties = rValues.GetMaterialProperties();
    CheckMaterialProperties(r_properties);
    const auto& r_options = rValues.GetOptions();
    auto& r_strain = rValues.GetStrainVector();
    if (r_options.IsNot(USE_ELEMENT_PROVIDED_STRAIN)) {
        const auto& r_f = rValues.GetDeformationGradientF();
        KRATOS_ERROR_IF(r_f.size1() != 3 || r_f.size2() != 3) << "Swift J2 requires a 3x3 F";
        r_strain.resize(6, false);
        for (SizeType i = 0; i < 3; ++i) {
            r_strain[i] = r_f(i, i) - 1.0;
        }
        r_strain[3] = r_f(0, 1) + r_f(1, 0);
        r_strain[4] = r_f(1, 2) + r_f(2, 1);
        r_strain[5] = r_f(0, 2) + r_f(2, 0);
    }
    KRATOS_ERROR_IF(r_strain.size() != 6) << "Swift J2 requires six engineering strain components";
    Vector strain = r_strain;
    AddInitialStrainVectorContribution(strain);
    for (const double value : strain) {
        KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Swift J2: non-finite strain";
    }
    const double young = r_properties[YOUNG_MODULUS];
    const double nu = r_properties[POISSON_RATIO];
    const double two_g = young / (1.0 + nu);
    const double bulk = young / (3.0 * (1.0 - 2.0 * nu));
    Matrix elastic_matrix = ZeroMatrix(6, 6);
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(elastic_matrix, young, nu);
    Vector stress = prod(elastic_matrix, strain - mPlasticStrain);
    AddInitialStressVectorContribution(stress);
    const double mean_stress = stress[0] / 3.0 + stress[1] / 3.0 + stress[2] / 3.0;
    Vector deviator = stress;
    for (SizeType i = 0; i < 3; ++i) {
        deviator[i] -= mean_stress;
    }
    double trial_norm = 0.0;
    for (SizeType i = 0; i < 6; ++i) {
        trial_norm = std::hypot(trial_norm, deviator[i] * (i < 3 ? 1.0 : std::sqrt(2.0)));
    }
    KRATOS_ERROR_IF(!std::isfinite(trial_norm) || !std::isfinite(mean_stress))
        << "Swift J2: non-finite trial stress";

    rPlasticStrain = mPlasticStrain;
    rP = mAccumulatedPlasticStrain;
    double yield_stress = 0.0;
    double hardening = 0.0;
    EvaluateHardening(r_properties, rP, yield_stress, hardening);
    const double c = std::sqrt(2.0 / 3.0);
    const double trial_residual = trial_norm - c * yield_stress;
    Matrix tangent = elastic_matrix;
    if (trial_residual > 0.0) {
        // R(gamma) = r_trial - 2G*gamma - c*sigma_y(p_n+c*gamma).
        // R' = -2G-c^2*H < 0. The elastic correction brackets the unique root.
        double lower = 0.0;
        double upper = trial_residual / two_g;
        double delta_gamma = 0.0;
        double residual = trial_residual;
        constexpr SizeType max_iterations = 100;
        const double stress_tolerance = 1.0e-12 * trial_norm;
        bool converged = false;
        for (SizeType iteration = 0; iteration < max_iterations; ++iteration) {
            rP = mAccumulatedPlasticStrain + c * delta_gamma;
            EvaluateHardening(r_properties, rP, yield_stress, hardening);
            residual = trial_norm - two_g * delta_gamma - c * yield_stress;
            if (std::abs(residual) <= stress_tolerance) {
                converged = true;
                break;
            }
            if (residual > 0.0) {
                lower = delta_gamma;
            } else {
                upper = delta_gamma;
            }
            const double derivative = two_g + c * c * hardening;
            const double candidate = delta_gamma + residual / derivative;
            delta_gamma = (candidate > lower && candidate < upper)
                ? candidate : lower + 0.5 * (upper - lower);
        }
        KRATOS_ERROR_IF_NOT(converged) << "Swift J2 radial return failed after " << max_iterations
            << " iterations: residual = " << residual << ", tolerance = " << stress_tolerance
            << ", committed p = " << mAccumulatedPlasticStrain << ", trial norm = " << trial_norm;

        // trial_norm is strictly positive here because the Swift yield stress is positive.
        const Vector normal = deviator / trial_norm;
        const double correction = two_g * delta_gamma / trial_norm;
        const double a = 1.0 - correction;
        for (SizeType i = 0; i < 6; ++i) {
            stress[i] = a * deviator[i] + (i < 3 ? mean_stress : 0.0);
            rPlasticStrain[i] += delta_gamma * normal[i] * (i < 3 ? 1.0 : 2.0);
        }
        if (r_options.Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
            // dr_trial = 2G*N:d_epsilon; dgamma = dr_trial/(2G+c^2*H).
            // ds = 2G*a*Pdev:d_epsilon - b*N*(N:d_epsilon),
            // b = 2G*(2G/(2G+c^2*H) - 2G*gamma/r_trial).
            // In engineering Voigt notation Pdev has shear diagonal 1/2;
            // N:d_epsilon is the ordinary dot product with engineering strain.
            // Consequently N_i*N_j needs no extra shear factors in the matrix.
            const double b = two_g * (1.0 / (1.0 + c * c * hardening / two_g) - correction);
            for (SizeType i = 0; i < 6; ++i) {
                for (SizeType j = 0; j < 6; ++j) {
                    const double volumetric = (i < 3 && j < 3) ? bulk : 0.0;
                    tangent(i, j) = volumetric + a * (elastic_matrix(i, j) - volumetric)
                        - b * normal[i] * normal[j];
                }
            }
        }
    }
    if (r_options.Is(COMPUTE_STRESS)) {
        rValues.GetStressVector() = stress;
    }
    if (r_options.Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
        rValues.GetConstitutiveMatrix() = tangent;
    }
    KRATOS_CATCH("")
}

bool SmallStrainSwiftJ2Plasticity3D::Has(const Variable<double>& rVariable)
{
    return rVariable == ACCUMULATED_PLASTIC_STRAIN;
}

bool SmallStrainSwiftJ2Plasticity3D::Has(const Variable<Vector>& rVariable)
{
    return rVariable == PLASTIC_STRAIN_VECTOR;
}

double& SmallStrainSwiftJ2Plasticity3D::GetValue(const Variable<double>& rVariable, double& rValue)
{
    if (rVariable == ACCUMULATED_PLASTIC_STRAIN) {
        rValue = mAccumulatedPlasticStrain;
        return rValue;
    }
    return BaseType::GetValue(rVariable, rValue);
}

Vector& SmallStrainSwiftJ2Plasticity3D::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
        return rValue;
    }
    return BaseType::GetValue(rVariable, rValue);
}

void SmallStrainSwiftJ2Plasticity3D::SetValue(const Variable<double>& rVariable,
    const double& rValue, const ProcessInfo& rProcessInfo)
{
    if (rVariable == ACCUMULATED_PLASTIC_STRAIN) {
        KRATOS_ERROR_IF(!std::isfinite(rValue) || rValue < 0.0) << "Swift J2: p must be finite and non-negative";
        mAccumulatedPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rVariable, rValue, rProcessInfo);
    }
}

void SmallStrainSwiftJ2Plasticity3D::SetValue(const Variable<Vector>& rVariable,
    const Vector& rValue, const ProcessInfo& rProcessInfo)
{
    if (rVariable == PLASTIC_STRAIN_VECTOR) {
        KRATOS_ERROR_IF(rValue.size() != 6) << "Swift J2: plastic strain must have six components";
        for (const double value : rValue) {
            KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Swift J2: plastic strain must be finite";
        }
        mPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rVariable, rValue, rProcessInfo);
    }
}

double& SmallStrainSwiftJ2Plasticity3D::CalculateValue(Parameters& rValues,
    const Variable<double>& rVariable, double& rValue)
{
    if (rVariable == YIELD_STRESS) {
        CheckMaterialProperties(rValues.GetMaterialProperties());
        double hardening = 0.0;
        EvaluateHardening(rValues.GetMaterialProperties(), mAccumulatedPlasticStrain, rValue, hardening);
        return rValue;
    }
    return GetValue(rVariable, rValue);
}

void SmallStrainSwiftJ2Plasticity3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("PlasticStrain", mPlasticStrain);
    rSerializer.save("AccumulatedPlasticStrain", mAccumulatedPlasticStrain);
}

void SmallStrainSwiftJ2Plasticity3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("PlasticStrain", mPlasticStrain);
    rSerializer.load("AccumulatedPlasticStrain", mAccumulatedPlasticStrain);
}

} // namespace Kratos
