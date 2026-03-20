//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Siemer
//

// System includes
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/beam_section_properties_utility.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesRod, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("ROD", std::vector<double>{0.04});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.005026548245743669, 1.0e-12);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.0106192982974676e-06, 1.0e-16);
    KRATOS_EXPECT_NEAR(section_properties.I33, 2.0106192982974676e-06, 1.0e-16);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 4.021238596594935e-06, 1.0e-16);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.9, 1.0e-12);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.9, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesTube, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("TUBE", std::vector<double>{0.05, 0.03});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.005026548245743671, 1.0e-12);
    KRATOS_EXPECT_NEAR(section_properties.I22, 4.2725660088821195e-06, 1.0e-16);
    KRATOS_EXPECT_NEAR(section_properties.I33, 4.2725660088821195e-06, 1.0e-16);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 8.545132017764239e-06, 1.0e-16);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.5, 1.0e-12);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.5, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesBar, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("BAR", std::vector<double>{0.05, 0.2});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.010000000000000002, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 3.333333333333334e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 2.083333333333334e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 0.0069653333333333355, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.8333333333333334, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.8333333333333334, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesBox, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("BOX", std::vector<double>{0.16, 0.24, 0.02, 0.015});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.012399999999999998, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 9.765333333333332e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 4.530333333333335e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 9.286174904942968e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.44354838709677424, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.38978494623655924, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesIProfileCaseInsensitive, KratosCoreFastSuite)
{
    const auto section_properties = BeamSectionPropertiesUtility::CalculateProperties(
        "i",
        std::vector<double>{0.20, 0.05, 0.05, 0.005, 0.005, 0.005});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0014500000000000001, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 7.6120833333333345e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 1.0614583333333336e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 1.2291666666666668e-08, 1.0e-20);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.5459770114942528, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.28735632183908044, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesChan, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("CHAN", std::vector<double>{0.08, 0.2, 0.01, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0048000000000000004, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.9440000000000003e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 3.0266666666666668e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 4.6000000000000004e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.2777777777777778, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.5555555555555556, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesT, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("T", std::vector<double>{0.12, 0.18, 0.03, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0066, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 1.9149545454545445e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 4.419999999999999e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 1.5199999999999998e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.3787878787878788, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.45454545454545453, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesCross, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("CROSS", std::vector<double>{0.18, 0.12, 0.02, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.006, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.0000000000000004e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.I33, 4.4999999999999996e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 0.00080048, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.5, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.3333333333333333, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesH, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("H", std::vector<double>{0.18, 0.12, 0.02, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.006, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.0000000000000004e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.I33, 4.4999999999999996e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 3.3599999999999996e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.5, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.3333333333333333, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesT1, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("T1", std::vector<double>{0.04, 0.08, 0.12, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0063999999999999994, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 6.933333333333335e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.I33, 1.8613333333333337e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 3.6266666666666674e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.625, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.20833333333333337, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesI1, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("I1", std::vector<double>{0.03, 0.02, 0.01, 0.08});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0037, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.1308333333333334e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 7.358333333333336e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 7.040000000000001e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.36036036036036034, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.472972972972973, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesChan1, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("CHAN1", std::vector<double>{0.03, 0.02, 0.01, 0.08});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0049, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 3.5370578231292517e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 9.743027210884355e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 8.650833333333332e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.272108843537415, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.5612244897959184, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesZ, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("Z", std::vector<double>{0.03, 0.02, 0.01, 0.08});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0037000000000000006, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.1308333333333334e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 1.5233333333333337e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 1.0708333333333336e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.36036036036036034, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.472972972972973, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesChan2, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("CHAN2", std::vector<double>{0.02, 0.03, 0.08, 0.06});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0038, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 2.0674561403508773e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 1.4066666666666662e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 4.859999999999999e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.4385964912280702, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.39473684210526316, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesT2, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("T2", std::vector<double>{0.1, 0.16, 0.03, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.005600000000000001, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 1.2800952380952384e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 2.586666666666667e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 1.2466666666666667e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.38690476190476186, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.44642857142857134, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesBox1, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("BOX1", std::vector<double>{0.08, 0.12, 0.02, 0.02, 0.02, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0063999999999999994, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 9.813333333333333e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 4.693333333333334e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 9e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.4166666666666667, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.41666666666666674, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesHexa, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("HEXA", std::vector<double>{0.02, 0.1, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.0016, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 4.666666666666668e-08, 1.0e-20);
    KRATOS_EXPECT_NEAR(section_properties.I33, 9.06666666666667e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 1.1540000000000003e-05, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.8333333333333334, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.8333333333333334, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesHat, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("HAT", std::vector<double>{0.08, 0.02, 0.12, 0.03});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.006, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 4.296000000000001e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 1.5799999999999994e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 6.933333333333335e-07, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.4444444444444445, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.38888888888888884, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityCalculatesHat1, KratosCoreFastSuite)
{
    const auto section_properties =
        BeamSectionPropertiesUtility::CalculateProperties("HAT1", std::vector<double>{0.14, 0.08, 0.06, 0.02, 0.02});

    KRATOS_EXPECT_NEAR(section_properties.Area, 0.009600000000000001, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.I22, 7.3533333333333355e-06, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.I33, 1.0560000000000002e-05, 1.0e-18);
    KRATOS_EXPECT_NEAR(section_properties.TorsionalInertia, 4.010666666666667e-06, 1.0e-19);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorY, 0.13888888888888887, 1.0e-15);
    KRATOS_EXPECT_NEAR(section_properties.ShearFactorZ, 0.6944444444444444, 1.0e-15);
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityThrowsOnUnsupportedSection, KratosCoreFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        BeamSectionPropertiesUtility::CalculateProperties("DBOX", std::vector<double>{1.0, 2.0}),
        "Unsupported CROSS_SECTION_TYPE");
}

KRATOS_TEST_CASE_IN_SUITE(BeamSectionPropertiesUtilityThrowsOnMissingDimension, KratosCoreFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        BeamSectionPropertiesUtility::CalculateProperties("BOX", std::vector<double>{0.16, 0.24, 0.02}),
        "Missing dimension 4");
}

} // namespace Kratos::Testing
