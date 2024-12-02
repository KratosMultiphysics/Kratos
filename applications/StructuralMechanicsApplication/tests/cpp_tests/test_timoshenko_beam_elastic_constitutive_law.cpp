// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "structural_mechanics_fast_suite.h"
#include "custom_constitutive/timoshenko_beam_elastic_constitutive_law.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimoshenkoBeamElasticConstitutiveLaw_CalculatesCauchyStressVector, KratosStructuralMechanicsFastSuite) {
    TimoshenkoBeamElasticConstitutiveLaw law;
    ConstitutiveLaw::Parameters parameters;
    constexpr auto axial_strain = 2.0e-3;
    constexpr auto curvature = 3.0e-2;
    constexpr auto shear_strain = 4.0e-3;
    constexpr auto number_of_strain_components = 3;
    auto strain_vector = Vector{number_of_strain_components};
    strain_vector[0] = axial_strain;
    strain_vector[1] = curvature;
    strain_vector[2] = shear_strain;
    parameters.SetStrainVector(strain_vector);
    auto stress_vector = Vector{ZeroVector{number_of_strain_components}};
    parameters.SetStressVector(stress_vector);

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    constexpr auto poissons_ratio = 0.2;
    properties[POISSON_RATIO] = poissons_ratio;
    constexpr auto area = 0.2;
    properties[CROSS_AREA] = area;
    constexpr auto moment_of_inertia = 0.5;
    properties[I33] = moment_of_inertia;
    constexpr auto effective_shear_area = 0.3;
    properties[AREA_EFFECTIVE_Y] = effective_shear_area;
    parameters.SetMaterialProperties(properties);
    const auto shear_modulus = youngs_modulus / (2.0 * (1.0 + poissons_ratio));

    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    // We'd like to actually call `CalculateMaterialResponsePK2`, but since that requires Gennady's
    // changes, we'll use `CalculateMaterialResponseCauchy` for now
    law.CalculateMaterialResponseCauchy(parameters);

    auto expected_stress_vector = Vector{number_of_strain_components};
    expected_stress_vector[0] = youngs_modulus * area * axial_strain;
    expected_stress_vector[1] = youngs_modulus * moment_of_inertia * curvature;
    expected_stress_vector[2] = shear_modulus * effective_shear_area * shear_strain;

    constexpr auto relative_tolerance = 1.0e-6;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(parameters.GetStressVector(), expected_stress_vector, relative_tolerance);
}

}