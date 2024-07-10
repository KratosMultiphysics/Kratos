// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
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
    Vector temp_strain = ScalarVector(1, 0.005);
    Vector temp_stress = ZeroVector(1);
    parameters.SetStrainVector(temp_strain);
    parameters.SetStressVector(temp_stress);

    Properties properties;
    properties[YOUNG_MODULUS] = 1e6;
    parameters.SetMaterialProperties(properties);
    law.CalculateMaterialResponsePK2(parameters);

    const double expected_stress = 5000; // = Strain * Young's Modulus
    KRATOS_EXPECT_EQ(expected_stress, parameters.GetStressVector()[0]);
}

KRATOS_TEST_CASE_IN_SUITE(TrussConstitutiveLaw_CalculatesLinearElasticStress_WithInitialStress, KratosStructuralMechanicsFastSuite) {
    TrussConstitutiveLaw law;
    ConstitutiveLaw::Parameters parameters;
    Vector temp_strain = ScalarVector(1, 0.005);
    Vector temp_stress = ZeroVector(1);
    parameters.SetStrainVector(temp_strain);
    parameters.SetStressVector(temp_stress);

    Properties properties;
    properties[YOUNG_MODULUS] = 1e6;
    parameters.SetMaterialProperties(properties);

    Vector initial_strain = ZeroVector(1);
    Vector initial_stress = ScalarVector(1, 5000);
    InitialState::Pointer p_initial_state = Kratos::make_intrusive<InitialState>(initial_strain, initial_stress);
    law.SetInitialState(p_initial_state);

    law.CalculateMaterialResponsePK2(parameters);

    constexpr double expected_stress = 10000; // = Strain * Young's Modulus + Initial Stress
    KRATOS_EXPECT_EQ(expected_stress, parameters.GetStressVector()[0]);
}

}
