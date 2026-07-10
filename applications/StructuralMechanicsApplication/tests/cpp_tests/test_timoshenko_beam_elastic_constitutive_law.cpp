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

namespace
{
using namespace Kratos;

constexpr auto number_of_strain_components = 3;

Vector MakeStrainVectorForTesting() {
    auto result = Vector{number_of_strain_components};
    result[0] = 2.0e-3; // axial strain
    result[1] = 3.0e-2; // curvature
    result[2] = 4.0e-3; // shear strain
    return result;
}

Properties MakePropertiesForTesting() {
    auto result = Properties{};
    result[YOUNG_MODULUS] = 1.0e6;
    result[POISSON_RATIO] = 0.2;
    result[CROSS_AREA] = 0.2;
    result[I33] = 0.5;
    result[AREA_EFFECTIVE_Y] = 0.3;
    return result;
}

Vector MakeInitialStrainVectorForTesting() {
    auto result = Vector{number_of_strain_components};
    result[0] = 4.0e-3; // initial axial strain
    result[1] = 5.0e-2; // initial curvature
    result[2] = 2.0e-3; // initial shear strain
    return result;
}

Vector MakeInitialStressVectorForTesting() {
    auto result = Vector{number_of_strain_components};
    result[0] = 5000.0;
    result[1] = 3000.0;
    result[2] = 2000.0;
    return result;
}

}

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimoshenkoBeamElasticConstitutiveLaw_CalculatesCauchyStressVector, KratosStructuralMechanicsFastSuite) {
    TimoshenkoBeamElasticConstitutiveLaw law;
    ConstitutiveLaw::Parameters parameters;
    auto strain_vector = MakeStrainVectorForTesting();
    parameters.SetStrainVector(strain_vector);
    auto stress_vector = Vector{ZeroVector{number_of_strain_components}};
    parameters.SetStressVector(stress_vector);

    const auto properties = MakePropertiesForTesting();
    parameters.SetMaterialProperties(properties);

    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    law.CalculateMaterialResponseCauchy(parameters);

    auto expected_stress_vector = Vector{number_of_strain_components};
    expected_stress_vector[0] = 400.0;
    expected_stress_vector[1] = 15000.0;
    expected_stress_vector[2] = 500.0;

    constexpr auto relative_tolerance = 1.0e-6;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(parameters.GetStressVector(), expected_stress_vector, relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(TimoshenkoBeamElasticConstitutiveLaw_CalculatesCauchyStressVector_WithInitialState, KratosStructuralMechanicsFastSuite) {
    TimoshenkoBeamElasticConstitutiveLaw law;
    ConstitutiveLaw::Parameters parameters;
    auto strain_vector = MakeStrainVectorForTesting();
    parameters.SetStrainVector(strain_vector);
    auto stress_vector = Vector{ZeroVector{number_of_strain_components}};
    parameters.SetStressVector(stress_vector);

    const auto properties = MakePropertiesForTesting();
    parameters.SetMaterialProperties(properties);

    const auto initial_strain_vector = MakeInitialStrainVectorForTesting();
    const auto initial_stress_vector = MakeInitialStressVectorForTesting();
    const auto p_initial_state = make_intrusive<InitialState>(
        initial_strain_vector, initial_stress_vector);
    law.SetInitialState(p_initial_state);

    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    law.CalculateMaterialResponseCauchy(parameters);

    auto expected_stress_vector = Vector{number_of_strain_components};
    expected_stress_vector[0] = 4600.0;
    expected_stress_vector[1] = -7000.0;
    expected_stress_vector[2] = 2250.0;

    constexpr auto relative_tolerance = 1.0e-6;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(parameters.GetStressVector(), expected_stress_vector, relative_tolerance);
}

}