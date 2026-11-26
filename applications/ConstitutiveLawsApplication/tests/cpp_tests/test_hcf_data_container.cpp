// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_constitutive/auxiliary_files/hcf_data_container.h"
#include "constitutive_laws_application_variables.h"

// Application includes
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"

namespace Kratos::Testing
{

/**
* Check the fatigue parameters
*/
KRATOS_TEST_CASE_IN_SUITE(HCFDataContainerFatigueTest, KratosConstitutiveLawsFastSuite)
{
    Properties material_properties;
    ProcessInfo process_info;
    HCFDataContainer mFatigueData = HCFDataContainer();
    HCFDataContainer::FatigueVariables HCFVariables = HCFDataContainer::FatigueVariables();

    Vector hcf_coefficients = ZeroVector(7);
    hcf_coefficients[0] = 0.322;
    hcf_coefficients[1] = 1.000;
    hcf_coefficients[2] = 0.6;
    hcf_coefficients[3] = 0.0309;
    hcf_coefficients[4] = 2.28815;
    hcf_coefficients[5] = -0.0284;
    hcf_coefficients[6] = 0.002;

    material_properties.SetValue(YOUNG_MODULUS, 201000000000);
    material_properties.SetValue(POISSON_RATIO, 0.263);
    material_properties.SetValue(FRACTURE_ENERGY, 4.1980284500E+06);
    material_properties.SetValue(YIELD_STRESS, 806.715E+06);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    material_properties.SetValue(HIGH_CYCLE_FATIGUE_COEFFICIENTS, hcf_coefficients);

    HCFVariables.MaxStress = 4.8419e+08;
    HCFVariables.MinStress = -4.84212e+08;
    HCFVariables.GlobalNumberOfCycles = 12;
    HCFVariables.LocalNumberOfCycles = 12;
    HCFVariables.ReversionFactor = -1.00005;
    HCFVariables.UltimateStress = 806.715E+06;

    mFatigueData.CalculateFatigueParameters(material_properties, HCFVariables);

    KRATOS_EXPECT_NEAR(HCFVariables.Alphat, 0.0309, 1.0e-6);
    KRATOS_EXPECT_NEAR(HCFVariables.Sth, 2.6071e+08, 1.0e+04);
    KRATOS_EXPECT_NEAR(HCFVariables.CyclesToFailure, 22406, 1.0);
    KRATOS_EXPECT_NEAR(HCFVariables.B0, 0.000231689, 1.0e-9);

    mFatigueData.CalculateFatigueReductionFactorAndWohlerStress(material_properties, HCFVariables);

    KRATOS_EXPECT_NEAR(HCFVariables.FatigueReductionFactor, 0.999655, 1.0e-6);
    KRATOS_EXPECT_NEAR(HCFVariables.WohlerStress, 0.975555, 1.0e-6);
}
} // namespace Kratos::Testing
