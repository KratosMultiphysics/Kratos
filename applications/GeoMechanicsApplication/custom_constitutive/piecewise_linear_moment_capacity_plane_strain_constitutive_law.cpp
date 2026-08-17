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
#include <limits>

using namespace std::string_literals;

namespace Kratos
{

ConstitutiveLaw::Pointer PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::Clone() const
{
    return make_shared<PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw>(*this);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize     = strain_size;
    rFeatures.mSpaceDimension = space_dimension;
}

bool PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::RequiresFinalizeMaterialResponse()
{
    return true;
}

double& PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateValue(Parameters& rParameters,
                                                                                const Variable<double>& rVariable,
                                                                                double& rValue)
{
    if (rVariable == TANGENT_MODULUS) {
        const auto curvature = rParameters.GetStrainVector()[1];
        rValue               = CalculateMomentAndTangentModulus(curvature).second;
        return rValue;
    }

    return BaseType::CalculateValue(rParameters, rVariable, rValue);
}

ConstitutiveLaw::StressMeasure PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::GetStressMeasure()
{
    return StressMeasure_PK2;
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rParameters)
{
    CalculateMaterialResponseCauchy(rParameters);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateMaterialResponseCauchy(Parameters& rParameters)
{
    const auto& r_options             = rParameters.GetOptions();
    auto&       r_material_properties = rParameters.GetMaterialProperties();
    auto&       r_strain_vector       = rParameters.GetStrainVector();

    const auto axial_strain = r_strain_vector[0]; // eps
    const auto curvature    = r_strain_vector[1]; // Kappa
    const auto shear_strain = r_strain_vector[2]; // Gamma_xy

    const auto   E  = r_material_properties[YOUNG_MODULUS];
    const auto   A  = r_material_properties[THICKNESS];
    const double I  = std::pow(A, 3) / 12.0; // Per unit length
    const auto   nu = r_material_properties[POISSON_RATIO];

    const auto G   = E / (2.0 * (1.0 + nu));
    const auto A_s = r_material_properties[THICKNESS_EFFECTIVE_Y];

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto& r_generalized_stress_vector = rParameters.GetStressVector();
        r_generalized_stress_vector.resize(strain_size, false);

        const auto   one_minus_nu_squared = 1.0 - nu * nu;
        const auto   EA_nu                = E * A / one_minus_nu_squared;
        const double EI_nu                = E * I / one_minus_nu_squared;
        const auto   GAs                  = G * A_s;

        const auto [moment, tangent_modulus] = CalculateMomentAndTangentModulus(curvature);

        r_generalized_stress_vector[0] = EA_nu * axial_strain; // Nx
        r_generalized_stress_vector[1] = moment;               // Mz
        r_generalized_stress_vector[2] = GAs * shear_strain;   // Vxy
        // Plane strain generalized stress components
        r_generalized_stress_vector[3] = nu * r_generalized_stress_vector[0]; // Nz
        r_generalized_stress_vector[4] = nu * r_generalized_stress_vector[1]; // Mx

        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto& r_stress_derivatives = rParameters.GetConstitutiveMatrix();
            r_stress_derivatives.resize(strain_size, strain_size, false);

            noalias(r_stress_derivatives) = ZeroMatrix(strain_size, strain_size);
            r_stress_derivatives(0, 0)    = EA_nu;           // dN_dEl
            r_stress_derivatives(1, 1)    = tangent_modulus; // dM_dkappa
            r_stress_derivatives(2, 2)    = GAs;             // dV_dGamma_xy
            r_stress_derivatives(3, 3)    = nu * EA_nu;      // dNz_dEl
            r_stress_derivatives(4, 4)    = nu * EI_nu;      // dMx_dkappa
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

    check_properties.CheckAvailabilityAndNotEmpty(GEO_PIECEWISE_LINEAR_MOMENT_LAW);

    // First provided point must be non-zero since (0,0) is implicitly added
    const auto& r_kappa_moments = rMaterialProperties[GEO_PIECEWISE_LINEAR_MOMENT_LAW];
    KRATOS_ERROR_IF(r_kappa_moments.size() == 0)
        << "GEO_PIECEWISE_LINEAR_MOMENT_LAW is empty for material with ID: " << rMaterialProperties.Id()
        << std::endl;
    KRATOS_ERROR_IF_NOT(r_kappa_moments[0].size() == 2)
        << "GEO_PIECEWISE_LINEAR_MOMENT_LAW has to have two components for material with ID: "
        << rMaterialProperties.Id() << std::endl;

    const auto& first_point  = r_kappa_moments[0];
    const auto  first_kappa  = first_point[0];
    const auto  first_moment = first_point[1];

    constexpr auto tolerance = std::numeric_limits<double>::min();
    KRATOS_ERROR_IF(std::abs(first_kappa) < tolerance || std::abs(first_moment) < tolerance)
        << "First provided point must be non-zero when assuming implicit (0,0). Got ("
        << first_kappa << ", " << first_moment << ")" << std::endl;

    Vector kappas(r_kappa_moments.size());
    Vector moments(r_kappa_moments.size());
    std::transform(r_kappa_moments.begin(), r_kappa_moments.end(), kappas.begin(),
                   [](const Vector& point) { return point[0]; });
    std::transform(r_kappa_moments.begin(), r_kappa_moments.end(), moments.begin(),
                   [](const Vector& point) { return point[1]; });
    CheckUtilities::CheckValuesAreAscending(kappas, "KAPPA in GEO_PIECEWISE_LINEAR_LAW");
    constexpr auto allow_equal = true;
    CheckUtilities::CheckValuesAreAscending(moments, "MOMENT in GEO_PIECEWISE_LINEAR_LAW", allow_equal);

    check_properties.Check(YOUNG_MODULUS);
    check_properties.Check(THICKNESS);
    check_properties.Check(THICKNESS_EFFECTIVE_Y);
    check_properties.Check(POISSON_RATIO, -1.0, 0.5);
    if (rMaterialProperties.Has(GEO_UNLOADING_RELOADING_MODULUS))
        check_properties.Check(GEO_UNLOADING_RELOADING_MODULUS);

    const auto nu = rMaterialProperties[POISSON_RATIO];
    const auto elastic_EI =
        rMaterialProperties[YOUNG_MODULUS] * rMaterialProperties[THICKNESS] / (1.0 - nu * nu);
    for (auto i = std::size_t{1}; i < kappas.size(); ++i) {
        const auto delta_kappa  = kappas[i] - kappas[i - 1];
        const auto delta_moment = moments[i] - moments[i - 1];
        const auto dM_dKappa    = delta_moment / delta_kappa;
        KRATOS_ERROR_IF(dM_dKappa >= elastic_EI)
            << "Derivative dM/dKappa must be smaller than elastic EI. Segment " << i
            << " has dM/dKappa = " << dM_dKappa << ", elastic EI = " << elastic_EI << std::endl;
        if (rMaterialProperties.Has(GEO_UNLOADING_RELOADING_MODULUS)) {
            KRATOS_ERROR_IF(dM_dKappa >= rMaterialProperties[GEO_UNLOADING_RELOADING_MODULUS])
                << "Derivative dM/dKappa must be smaller than unloading/reloading modulus. Segment "
                << i << " has dM/dKappa = " << dM_dKappa << ", unloading/reloading modulus = "
                << rMaterialProperties[GEO_UNLOADING_RELOADING_MODULUS] << std::endl;
        }
    }

    return 0;
}

std::string PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::Info() const
{
    return "PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw"s;
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("BendingMomentCurvatureTable", mBendingMomentCurvatureTable);
    rSerializer.save("AccumulatedCurvature", mAccumulatedCurvature);
    rSerializer.save("UnReloadCenter", mUnReLoadCenter);
    rSerializer.save("UnReloadModulus", mUnReLoadModulus);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("BendingMomentCurvatureTable", mBendingMomentCurvatureTable);
    rSerializer.load("AccumulatedCurvature", mAccumulatedCurvature);
    rSerializer.load("UnReloadCenter", mUnReLoadCenter);
    rSerializer.load("UnReloadModulus", mUnReLoadModulus);
}

void PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties, const ConstitutiveLaw::GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
{
    BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    mBendingMomentCurvatureTable.Clear();

    const auto& r_kappa_moments = rMaterialProperties[GEO_PIECEWISE_LINEAR_MOMENT_LAW];
    // include implicit origin (0,0) as the first table point
    mBendingMomentCurvatureTable.PushBack(0.0, 0.0);
    for (auto i = std::size_t{0}; i < r_kappa_moments.size(); ++i) {
        mBendingMomentCurvatureTable.PushBack(r_kappa_moments[i][0], r_kappa_moments[i][1]);
    }
    // Append a final plateau point beyond the last provided kappa so the
    // piecewise-linear table remains constant for larger curvatures
    const std::size_t last_index = r_kappa_moments.size() - 1;
    mBendingMomentCurvatureTable.PushBack(r_kappa_moments[last_index][0] + 1.0,
                                          r_kappa_moments[last_index][1]);

    // Reset unload/reload state
    mAccumulatedCurvature = 0.0;
    mUnReLoadCenter       = 0.0;
    // Store optional modulus if provided
    mUnReLoadModulus = rMaterialProperties.Has(GEO_UNLOADING_RELOADING_MODULUS)
                           ? rMaterialProperties[GEO_UNLOADING_RELOADING_MODULUS]
                           : 0.0;
}

std::pair<double, double> PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::CalculateMomentAndTangentModulus(double Curvature) const
{
    auto moment          = 0.0;
    auto tangent_modulus = 0.0;

    if (mUnReLoadModulus > 0.0) {
        if (IsWithinUnReLoading(Curvature)) {
            moment          = mUnReLoadModulus * (Curvature - mUnReLoadCenter);
            tangent_modulus = mUnReLoadModulus;
        } else {
            const auto curvature_amplitude = CalculateUnReLoadAmplitude();
            const auto effective_curvature =
                mAccumulatedCurvature + (std::abs(Curvature - mUnReLoadCenter) - curvature_amplitude);
            moment          = BackboneMoment(effective_curvature);
            moment          = Curvature - mUnReLoadCenter < 0.0 ? -moment : moment;
            tangent_modulus = BackboneTangentModulus(effective_curvature);
        }
    } else {
        moment          = BackboneMoment(Curvature);
        moment          = Curvature < 0.0 ? -moment : moment;
        tangent_modulus = BackboneTangentModulus(Curvature);
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

double PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::BackboneMoment(double Curvature) const
{
    return mBendingMomentCurvatureTable.GetValue(std::abs(Curvature));
}

double PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::BackboneTangentModulus(double Curvature) const
{
    return mBendingMomentCurvatureTable.GetDerivative(std::abs(Curvature));
}
} // namespace Kratos
