// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Gennady Markelov

// Project includes
#include "piecewise_linear_moment_capacity_constitutive_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"

namespace Kratos
{
using namespace std::string_literals;

ConstitutiveLaw::Pointer PiecewiseLinearMomentCapacityConstitutiveLaw::Clone() const
{
    return make_shared<PiecewiseLinearMomentCapacityConstitutiveLaw>(*this);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize     = strain_size;
    rFeatures.mSpaceDimension = space_dimenstion;
}

double& PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters,
                                                                     const Variable<double>& rVariable,
                                                                     double& rValue)
{
    if (rVariable == TANGENT_MODULUS) {
        const auto curvature = rParameters.GetStrainVector()[1];
        // If a modulus for un/reload is supplied, use the stored modulus inside the elastic window
        if (mUnReLoadModulus > 0.0) {
            // Use actual curvature (with sign) for unload/reload window checks
            if (IsWithinUnReLoading(curvature)) {
                rValue = mUnReLoadModulus;
                return rValue;
            }
            // outside the window: compute effective curvature on backbone (positive quantity)
            const auto amplitude = CalculateUnReLoadAmplitude();
            const auto effective = mAccumulatedCurvature + (std::abs(curvature - mUnReLoadCenter) - amplitude);
            rValue = mStressStrainTable.GetDerivative(effective);
            return rValue;
        }
        // Default: TANGENT_MODULUS is the derivative dM/d(kappa) based on absolute curvature
        rValue = mStressStrainTable.GetDerivative(std::abs(curvature));
        return rValue;
    }

    return BaseType::CalculateValue(rParameters, rVariable, rValue);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_law_options      = rValues.GetOptions();
    auto&       r_material_properties = rValues.GetMaterialProperties();
    auto&       r_strain_vector       = rValues.GetStrainVector();
    AddInitialStrainVectorContribution(r_strain_vector);

    const double axial_strain = r_strain_vector[0]; // E_l
    const double curvature    = r_strain_vector[1]; // Kappa
    const double shear_strain = r_strain_vector[2]; // Gamma_xy

    const double E  = r_material_properties[YOUNG_MODULUS];
    const double A  = r_material_properties[THICKNESS];
    const double nu = r_material_properties[POISSON_RATIO];

    const double G   = E / (2.0 * (1.0 + nu));
    const double A_s = r_material_properties[THICKNESS_EFFECTIVE_Y]; // Per unit length

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto& r_generalized_stress_vector = rValues.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);

        const double one_minus_nu_squared = 1.0 - nu * nu;
        const double EA_nu                = E * A / one_minus_nu_squared;
        const double GAs                  = G * A_s;

        r_generalized_stress_vector[0] = EA_nu * axial_strain; // Nx

        // Moment calculation using piecewise table with unload/reload logic
        double moment = 0.0;
        if (mUnReLoadModulus > 0.0) {
            if (IsWithinUnReLoading(curvature)) {
                moment = mUnReLoadModulus * (curvature - mUnReLoadCenter);
            } else {
                const auto amplitude = CalculateUnReLoadAmplitude();
                const auto effective =
                    mAccumulatedCurvature + (std::abs(curvature - mUnReLoadCenter) - amplitude);
                moment = mStressStrainTable.GetValue(effective);
                moment = curvature - mUnReLoadCenter < 0.0 ? -moment : moment;
            }
        } else {
            moment = mStressStrainTable.GetValue(std::abs(curvature));
            moment = curvature < 0.0 ? -moment : moment;
        }
        r_generalized_stress_vector[1] = moment;             // Mz
        r_generalized_stress_vector[2] = GAs * shear_strain; // Vxy

        AddInitialStressVectorContribution(r_generalized_stress_vector);

        if (r_material_properties.Has(BEAM_PRESTRESS_PK2)) {
            r_generalized_stress_vector += r_material_properties[BEAM_PRESTRESS_PK2];
        }

        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto& r_stress_derivatives = rValues.GetConstitutiveMatrix();
            if (r_stress_derivatives.size1() != strain_size || r_stress_derivatives.size2() != strain_size)
                r_stress_derivatives.resize(strain_size, strain_size, false);
            noalias(r_stress_derivatives) = ZeroMatrix(strain_size, strain_size);

            // Tangent modulus for moment
            double tangent_modulus = 0.0;
            if (mUnReLoadModulus > 0.0) {
                if (IsWithinUnReLoading(curvature)) {
                    tangent_modulus = mUnReLoadModulus;
                } else {
                    const auto amplitude = CalculateUnReLoadAmplitude();
                    const auto effective =
                        mAccumulatedCurvature + (std::abs(curvature - mUnReLoadCenter) - amplitude);
                    tangent_modulus = mStressStrainTable.GetDerivative(effective);
                }
            } else {
                tangent_modulus = mStressStrainTable.GetDerivative(std::abs(curvature));
            }

            r_stress_derivatives(0, 0) = EA_nu;           // dN_dEl
            r_stress_derivatives(1, 1) = tangent_modulus; // dM_dkappa
            r_stress_derivatives(2, 2) = GAs;             // dV_dGamma_xy
        }
    }
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    // compute stress first
    CalculateMaterialResponsePK2(rValues);

    const auto curvature = rValues.GetStrainVector()[1];

    if (mUnReLoadModulus > 0.0) {
        // Check for sign reversal: if sign of curvature changes, handle specially
        const bool sign_reversed =
            mPreviousCurvature != 0.0 && std::signbit(curvature) != std::signbit(mPreviousCurvature);

        if (sign_reversed) {
            // On sign reversal, reset center and accumulated to stay on the backbone
            // at the current absolute position (not overshooting)
            mAccumulatedCurvature = std::abs(curvature);
            const auto amplitude  = CalculateUnReLoadAmplitude();
            mUnReLoadCenter       = std::abs(curvature) > amplitude
                                        ? (curvature > 0.0 ? std::abs(curvature) - amplitude
                                                           : -(std::abs(curvature) - amplitude))
                                        : 0.0;
        } else if (!IsWithinUnReLoading(curvature)) {
            // Normal unload/reload update (no sign reversal)
            const auto difference = curvature - mUnReLoadCenter;
            const auto amplitude  = CalculateUnReLoadAmplitude();
            mAccumulatedCurvature += std::abs(difference) - amplitude;
            mUnReLoadCenter = difference > 0.0 ? curvature - amplitude : curvature + amplitude;
        }
    }
    mPreviousCurvature = curvature;
}

SizeType PiecewiseLinearMomentCapacityConstitutiveLaw::GetStrainSize() const { return strain_size; }

int PiecewiseLinearMomentCapacityConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const ProcessInfo& rCurrentProcessInfo) const
{
    if (const auto exit_code = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
        exit_code != 0)
        return exit_code;

    // Piecewise linear table properties
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(KAPPA_PIECEWISE_LINEAR_LAW))
        << "No KAPPA_PIECEWISE_LINEAR_LAW found" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MOMENTUM_PIECEWISE_LINEAR_LAW))
        << "No MOMENTUM_PIECEWISE_LINEAR_LAW found" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[KAPPA_PIECEWISE_LINEAR_LAW].empty())
        << "KAPPA_PIECEWISE_LINEAR_LAW is empty";
    KRATOS_ERROR_IF_NOT(rMaterialProperties[KAPPA_PIECEWISE_LINEAR_LAW].size() ==
                        rMaterialProperties[MOMENTUM_PIECEWISE_LINEAR_LAW].size())
        << "The number of strain components does not match the number of momentum components" << std::endl;

    // Since (0,0) is assumed, the first provided point must be non-zero
    const auto     first_kappa    = rMaterialProperties[KAPPA_PIECEWISE_LINEAR_LAW][0];
    const auto     first_momentum = rMaterialProperties[MOMENTUM_PIECEWISE_LINEAR_LAW][0];
    constexpr auto tolerance      = 1.0e-15;
    KRATOS_ERROR_IF(std::abs(first_kappa) < tolerance || std::abs(first_momentum) < tolerance)
        << "First provided point must be non-zero when assuming implicit (0,0). Got ("
        << first_kappa << ", " << first_momentum << ")" << std::endl;

    CheckUtilities::CheckValuesAreAscending(rMaterialProperties[KAPPA_PIECEWISE_LINEAR_LAW],
                                            "KAPPA_PIECEWISE_LINEAR_LAW");
    constexpr auto allow_equal = true;
    CheckUtilities::CheckValuesAreAscending(rMaterialProperties[MOMENTUM_PIECEWISE_LINEAR_LAW],
                                            "MOMENTUM_PIECEWISE_LINEAR_LAW", allow_equal);

    // Validate optional UNRELOAD_MODULUS if provided
    if (rMaterialProperties.Has(UNRELOAD_MODULUS)) {
        const auto value = rMaterialProperties[UNRELOAD_MODULUS];
        KRATOS_ERROR_IF(value <= 0.0) << "UNRELOAD_MODULUS must be positive. Got " << value << std::endl;
    }

    return 0;
}

std::string PiecewiseLinearMomentCapacityConstitutiveLaw::Info() const
{
    return "PiecewiseLinearMomentCapacityConstitutiveLaw"s;
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("StressStrainTable", mStressStrainTable);
    rSerializer.save("AccumulatedCurvature", mAccumulatedCurvature);
    rSerializer.save("PreviousCurvature", mPreviousCurvature);
    rSerializer.save("UnReloadCenter", mUnReLoadCenter);
    rSerializer.save("UnReloadModulus", mUnReLoadModulus);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("StressStrainTable", mStressStrainTable);
    rSerializer.load("AccumulatedCurvature", mAccumulatedCurvature);
    rSerializer.load("PreviousCurvature", mPreviousCurvature);
    rSerializer.load("UnReloadCenter", mUnReLoadCenter);
    rSerializer.load("UnReloadModulus", mUnReLoadModulus);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                                      const ConstitutiveLaw::GeometryType& rElementGeometry,
                                                                      const Vector& rShapeFunctionsValues)
{
    BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    mStressStrainTable.Clear();

    const auto& r_kappa   = rMaterialProperties[KAPPA_PIECEWISE_LINEAR_LAW];
    const auto& r_moments = rMaterialProperties[MOMENTUM_PIECEWISE_LINEAR_LAW];
    // include implicit origin (0,0) as the first table point
    mStressStrainTable.PushBack(0.0, 0.0);
    for (auto i = std::size_t{0}; i < r_kappa.size(); ++i) {
        mStressStrainTable.PushBack(r_kappa[i], r_moments[i]);
    }
    // Append a final plateau point beyond the last provided kappa so the
    // piecewise-linear table remains constant for larger curvatures
    const std::size_t last_index  = r_kappa.size() - 1;
    const auto        last_kappa  = r_kappa[last_index];
    const auto        last_moment = r_moments[last_index];
    mStressStrainTable.PushBack(last_kappa + 1.0, last_moment);

    // Reset unload/reload state
    mAccumulatedCurvature = 0.0;
    mPreviousCurvature    = 0.0;
    mUnReLoadCenter       = 0.0;
    // Store optional modulus if provided
    if (rMaterialProperties.Has(UNRELOAD_MODULUS)) {
        mUnReLoadModulus = rMaterialProperties[UNRELOAD_MODULUS];
    } else {
        mUnReLoadModulus = 0.0;
    }
}

double PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateUnReLoadAmplitude() const
{
    // amplitude is equal to M(accumulated) divided by modulus
    return mStressStrainTable.GetValue(mAccumulatedCurvature) / mUnReLoadModulus;
}

bool PiecewiseLinearMomentCapacityConstitutiveLaw::IsWithinUnReLoading(double Curvature) const
{
    return std::abs(Curvature - mUnReLoadCenter) < CalculateUnReLoadAmplitude();
}

} // namespace Kratos
