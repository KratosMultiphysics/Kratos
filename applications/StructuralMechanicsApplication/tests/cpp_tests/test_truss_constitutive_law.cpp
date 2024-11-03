// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Richard Faasse
//

#include "structural_mechanics_fast_suite.h"
#include "custom_constitutive/truss_constitutive_law.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TrussConstitutiveLaw_CalculatesLinearElasticStress, KratosStructuralMechanicsFastSuite) {
    TrussConstitutiveLaw law;
    ConstitutiveLaw::Parameters parameters;
    constexpr auto induced_strain = 0.005;
    Vector temp_strain = ScalarVector(1, induced_strain);
    Vector temp_stress = ZeroVector(1);
    parameters.SetStrainVector(temp_strain);
    parameters.SetStressVector(temp_stress);

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    parameters.SetMaterialProperties(properties);
    law.CalculateMaterialResponsePK2(parameters);

    constexpr double expected_stress = induced_strain * youngs_modulus;
    KRATOS_EXPECT_EQ(expected_stress, parameters.GetStressVector()[0]);
}

KRATOS_TEST_CASE_IN_SUITE(TrussConstitutiveLaw_CalculatesLinearElasticStress_WithInitialState, KratosStructuralMechanicsFastSuite) {
    TrussConstitutiveLaw law;
    constexpr auto induced_strain = 0.005;
    Vector temp_strain = ScalarVector(1, induced_strain);
    Vector temp_stress = ZeroVector(1);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(temp_strain);
    parameters.SetStressVector(temp_stress);

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    parameters.SetMaterialProperties(properties);

    constexpr auto initial_stress = 5000.0;
    constexpr auto initial_strain = 0.01;
    const auto p_initial_state = Kratos::make_intrusive<InitialState>(
        ScalarVector(1, initial_strain), ScalarVector(1, initial_stress));
    law.SetInitialState(p_initial_state);

    law.CalculateMaterialResponsePK2(parameters);

    constexpr double expected_stress =
        (induced_strain - initial_strain) * youngs_modulus + initial_stress;
    KRATOS_EXPECT_EQ(expected_stress, parameters.GetStressVector()[0]);
}

}
