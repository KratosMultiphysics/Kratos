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
#include "piecewise_linear_moment_capacity_plane_strain_constitutive_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"

namespace Kratos
{
using namespace std::string_literals;

ConstitutiveLaw::Pointer PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::Clone() const
{
    return make_shared<PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw>(*this);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize     = strain_size;
    rFeatures.mSpaceDimension = space_dimenstion;
}

bool PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::RequiresFinalizeMaterialResponse()
{
    return true;
}

double& PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameters, const Variable<double>& rVariable, double& rValue)
{
    if (rVariable == TANGENT_MODULUS) {
        const auto curvature                 = rParameters.GetStrainVector()[1];
        const auto [moment, tangent_modulus] = CalculateMomentAndTangentModulus(curvature);
        rValue                               = tangent_modulus;
        return rValue;
    }

    return BaseType::CalculateValue(rParameters, rVariable, rValue);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rParameters)
{
    CalculateMaterialResponseCauchy(rParameters);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_options             = rParameters.GetOptions();
    auto&       r_material_properties = rParameters.GetMaterialProperties();
    auto&       r_strain_vector       = rParameters.GetStrainVector();
    AddInitialStrainVectorContribution(r_strain_vector);

    const auto axial_strain = r_strain_vector[0]; // eps
    const auto curvature    = r_strain_vector[1]; // Kappa
    const auto shear_strain = r_strain_vector[2]; // Gamma_xy

    const auto E  = r_material_properties[YOUNG_MODULUS];
    const auto A  = r_material_properties[THICKNESS];
    const auto nu = r_material_properties[POISSON_RATIO];

    const auto G   = E / (2.0 * (1.0 + nu));
    const auto A_s = r_material_properties[THICKNESS_EFFECTIVE_Y];

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto& r_generalized_stress_vector = rParameters.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);

        const auto one_minus_nu_squared = 1.0 - nu * nu;
        const auto EA_nu                = E * A / one_minus_nu_squared;
        const auto GAs                  = G * A_s;

        const auto [moment, tangent_modulus] = CalculateMomentAndTangentModulus(curvature);

        r_generalized_stress_vector[0] = EA_nu * axial_strain; // Nx
        r_generalized_stress_vector[1] = moment;               // Mz
        r_generalized_stress_vector[2] = GAs * shear_strain;   // Vxy

        AddInitialStressVectorContribution(r_generalized_stress_vector);

        if (r_material_properties.Has(BEAM_PRESTRESS_PK2)) {
            r_generalized_stress_vector += r_material_properties[BEAM_PRESTRESS_PK2];
        }

        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto& r_stress_derivatives = rParameters.GetConstitutiveMatrix();
            if (r_stress_derivatives.size1() != strain_size || r_stress_derivatives.size2() != strain_size)
                r_stress_derivatives.resize(strain_size, strain_size, false);

            noalias(r_stress_derivatives) = ZeroMatrix(strain_size, strain_size);
            r_stress_derivatives(0, 0)    = EA_nu;           // dN_dEl
            r_stress_derivatives(1, 1)    = tangent_modulus; // dM_dkappa
            r_stress_derivatives(2, 2)    = GAs;             // dV_dGamma_xy
        }
    }
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rParameters)
{
    // compute stress first
    CalculateMaterialResponsePK2(rParameters);

    const auto curvature = rParameters.GetStrainVector()[1];

    if (mUnReLoadModulus > 0.0 && !IsWithinUnReLoading(curvature)) {
        const auto difference = curvature - mUnReLoadCenter;
        mAccumulatedCurvature += std::abs(difference) - CalculateUnReLoadAmplitude();
        mUnReLoadCenter = difference > 0.0 ? curvature - CalculateUnReLoadAmplitude()
                                           : curvature + CalculateUnReLoadAmplitude();
    }
}

SizeType PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::GetStrainSize() const
{
    return strain_size;
}

int PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::Check(const Properties& rMaterialProperties,
                                                                   const GeometryType& rElementGeometry,
                                                                   const ProcessInfo& rCurrentProcessInfo) const
{
    if (const auto exit_code = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
        exit_code != 0)
        return exit_code;

    CheckProperties check_properties(rMaterialProperties, "material properties",
                                     CheckProperties::Bounds::AllExclusive);

    check_properties.CheckAvailabilityAndNotEmpty(KAPPA_PIECEWISE_LINEAR_LAW);
    check_properties.CheckAvailabilityAndNotEmpty(MOMENTUM_PIECEWISE_LINEAR_LAW);

    const auto& r_kappa   = rMaterialProperties[KAPPA_PIECEWISE_LINEAR_LAW];
    const auto& r_moments = rMaterialProperties[MOMENTUM_PIECEWISE_LINEAR_LAW];
    KRATOS_ERROR_IF(r_kappa.size() != r_moments.size())
        << "The number of entries in KAPPA_PIECEWISE_LINEAR_LAW (" << r_kappa.size()
        << ") does not match MOMENTUM_PIECEWISE_LINEAR_LAW (" << r_moments.size() << ")" << std::endl;

    // First provided point must be non-zero since (0,0) is implicitly added
    constexpr auto tolerance      = 1.0e-15;
    const auto     first_kappa    = r_kappa[0];
    const auto     first_momentum = r_moments[0];
    KRATOS_ERROR_IF(std::abs(first_kappa) < tolerance || std::abs(first_momentum) < tolerance)
        << "First provided point must be non-zero when assuming implicit (0,0). Got ("
        << first_kappa << ", " << first_momentum << ")" << std::endl;

    CheckUtilities::CheckValuesAreAscending(r_kappa, "KAPPA_PIECEWISE_LINEAR_LAW");
    constexpr auto allow_equal = true;
    CheckUtilities::CheckValuesAreAscending(r_moments, "MOMENTUM_PIECEWISE_LINEAR_LAW", allow_equal);

    check_properties.Check(YOUNG_MODULUS);
    check_properties.Check(THICKNESS);
    check_properties.Check(THICKNESS_EFFECTIVE_Y);
    check_properties.Check(POISSON_RATIO, -1.0, 0.5);

    if (rMaterialProperties.Has(UNRELOAD_MODULUS)) check_properties.Check(UNRELOAD_MODULUS);

    return 0;
}

std::string PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::Info() const
{
    return "PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw"s;
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("StressStrainTable", mStressStrainTable);
    rSerializer.save("AccumulatedCurvature", mAccumulatedCurvature);
    rSerializer.save("UnReloadCenter", mUnReLoadCenter);
    rSerializer.save("UnReloadModulus", mUnReLoadModulus);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("StressStrainTable", mStressStrainTable);
    rSerializer.load("AccumulatedCurvature", mAccumulatedCurvature);
    rSerializer.load("UnReloadCenter", mUnReLoadCenter);
    rSerializer.load("UnReloadModulus", mUnReLoadModulus);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties, const ConstitutiveLaw::GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
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
    mUnReLoadCenter       = 0.0;
    // Store optional modulus if provided
    if (rMaterialProperties.Has(UNRELOAD_MODULUS)) {
        mUnReLoadModulus = rMaterialProperties[UNRELOAD_MODULUS];
    } else {
        mUnReLoadModulus = 0.0;
    }
}

std::pair<double, double> PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateMomentAndTangentModulus(double curvature) const
{
    double moment          = 0.0;
    double tangent_modulus = 0.0;

    if (mUnReLoadModulus > 0.0) {
        if (IsWithinUnReLoading(curvature)) {
            moment          = mUnReLoadModulus * (curvature - mUnReLoadCenter);
            tangent_modulus = mUnReLoadModulus;
        } else {
            const auto amplitude = CalculateUnReLoadAmplitude();
            const auto effective = mAccumulatedCurvature + (std::abs(curvature - mUnReLoadCenter) - amplitude);
            moment          = BackboneMoment(effective);
            moment          = curvature - mUnReLoadCenter < 0.0 ? -moment : moment;
            tangent_modulus = BackboneTangentModulus(effective);
        }
    } else {
        moment          = BackboneMoment(curvature);
        moment          = curvature < 0.0 ? -moment : moment;
        tangent_modulus = BackboneTangentModulus(curvature);
    }

    return {moment, tangent_modulus};
}

double PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateUnReLoadAmplitude() const
{
    return BackboneMoment(mAccumulatedCurvature) / mUnReLoadModulus;
}

bool PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::IsWithinUnReLoading(double Curvature) const
{
    return std::abs(Curvature - mUnReLoadCenter) < CalculateUnReLoadAmplitude();
}

double PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::BackboneMoment(double curvature) const
{
    return mStressStrainTable.GetValue(std::abs(curvature));
}

double PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::BackboneTangentModulus(double curvature) const
{
    return mStressStrainTable.GetDerivative(std::abs(curvature));
}
} // namespace Kratos
