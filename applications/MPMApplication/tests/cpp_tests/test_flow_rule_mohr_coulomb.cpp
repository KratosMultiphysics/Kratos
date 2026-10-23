//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes
#include <array>
#include <cmath>

// External includes

// Project includes
#include "testing/testing.h"

// Application includes
#include "mpm_application_variables.h"

// Hardening law
#include "custom_constitutive/hardening_laws/mpm_hardening_law.hpp"

// Yields
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"

// Flow Rule
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"


namespace Kratos
{
namespace Testing
{

    typedef Node NodeType;

    typedef MPMHardeningLaw HL;

    typedef MPMYieldCriterion YC;

    typedef MCYieldCriterion MCYC;

    typedef MPMFlowRule FR;

    typedef MCPlasticFlowRule MCFR;

    void GenerateTestMCVariables(
        Matrix& rStress, Matrix& rStrain,
        Properties &rMaterialProperties)
    {
        rStress = ZeroMatrix(3,3);
        rStrain = ZeroMatrix(3,3);
        rStrain(0,0) = 5.5e-2;
        rStrain(1,1) = -12.2e-2;
        rStrain(2,2) = 8.3e-2;

        rMaterialProperties.SetValue(DENSITY, 2200.0);
        rMaterialProperties.SetValue(YOUNG_MODULUS, 2000000.0);
        rMaterialProperties.SetValue(POISSON_RATIO, 0.3);
        rMaterialProperties.SetValue(INTERNAL_FRICTION_ANGLE, 0.5235987755982);
        rMaterialProperties.SetValue(COHESION, 5000.0);
        rMaterialProperties.SetValue(INTERNAL_DILATANCY_ANGLE, 0.436332312998582);

    }

    /**
    * Check the correct calculation of the integrated stress with the CL's
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulomb, KratosMPMFastSuite)
    {
        Matrix stress, strain;
        Properties material_properties;
        GenerateTestMCVariables(stress, strain, material_properties);

        // Construct Flow Rule
        HL::Pointer mc_hl_pointer = HL::Pointer( new HL() );
        YC::Pointer mc_yc_pointer = YC::Pointer( new MCYC() );
        FR::Pointer mc_fr_pointer = FR::Pointer( new MCFR() );
        mc_fr_pointer->InitializeMaterial(mc_yc_pointer, mc_hl_pointer, material_properties);

        // Compute trial elastic stresses
        FR::RadialReturnVariables rma_variables;
        mc_fr_pointer->CalculatePrincipalStressTrial(rma_variables, strain, stress, material_properties);

        Vector stress_trial_analytic = ZeroVector(3);
        stress_trial_analytic[0] =  1.03076923077e+05;
        stress_trial_analytic[1] = -1.69230769230e+05;
        stress_trial_analytic[2] =  1.46153846154e+05;

        KRATOS_EXPECT_NEAR(stress(0,0), stress_trial_analytic[0], 1e-6);
        KRATOS_EXPECT_NEAR(stress(1,1), stress_trial_analytic[1], 1e-6);
        KRATOS_EXPECT_NEAR(stress(2,2), stress_trial_analytic[2], 1e-6);

        // Compute new stresses after return mapping
        Matrix dummy_deformation_gradient = IdentityMatrix(3);
        mc_fr_pointer->CalculateReturnMapping( rma_variables, dummy_deformation_gradient, stress, strain, material_properties);

        Vector stress_analytic = ZeroVector(3);
        stress_analytic[0] = -5.31271630223e+04;
        stress_analytic[1] = -1.76701997143e+05;
        stress_analytic[2] = -5.31271630223e+04;
        const double yield_analytic = 5.90371799617e+05;

        KRATOS_EXPECT_NEAR(rma_variables.TrialStateFunction, yield_analytic, 1e-6);
        KRATOS_EXPECT_NEAR(stress(0,0), stress_analytic[0], 1e-6);
        KRATOS_EXPECT_NEAR(stress(1,1), stress_analytic[1], 1e-6);
        KRATOS_EXPECT_NEAR(stress(2,2), stress_analytic[2], 1e-6);

    }

    namespace
    {
        void check_mc_return_mapping(
            const std::array<double,3>& rTrialStress,
            const std::array<double,3>& rExpectedStress,
            const unsigned int expected_region,
            const double friction_angle = 0.0,
            const double dilatancy_angle = 0.0)
        {
            Matrix stress, strain;
            Properties material_properties;
            GenerateTestMCVariables(stress, strain, material_properties);
            material_properties.SetValue(INTERNAL_FRICTION_ANGLE, friction_angle);
            material_properties.SetValue(INTERNAL_DILATANCY_ANGLE, dilatancy_angle);

            // Invert isotropic elasticity to prescribe the trial principal stresses.
            const double young_modulus = material_properties[YOUNG_MODULUS];
            const double poisson_ratio = material_properties[POISSON_RATIO];
            const double trial_trace = rTrialStress[0] + rTrialStress[1] + rTrialStress[2];
            for (unsigned int i = 0; i < 3; ++i) {
                strain(i,i) = ((1.0 + poisson_ratio) * rTrialStress[i] - poisson_ratio * trial_trace) / young_modulus;
            }

            HL::Pointer p_hardening_law = Kratos::make_shared<HL>();
            YC::Pointer p_yield_criterion = Kratos::make_shared<MCYC>();
            FR::Pointer p_flow_rule = Kratos::make_shared<MCFR>();
            p_flow_rule->InitializeMaterial(p_yield_criterion, p_hardening_law, material_properties);

            FR::RadialReturnVariables return_mapping_variables;
            return_mapping_variables.initialize();
            p_flow_rule->CalculatePrincipalStressTrial(return_mapping_variables, strain, stress, material_properties);
            for (unsigned int i = 0; i < 3; ++i) {
                KRATOS_EXPECT_NEAR(stress(i,i), rTrialStress[i], 1e-8);
            }

            const Matrix deformation_gradient = IdentityMatrix(3);
            const bool plasticity_active = p_flow_rule->CalculateReturnMapping(
                return_mapping_variables, deformation_gradient, stress, strain, material_properties);
            KRATOS_EXPECT_EQ(plasticity_active, expected_region != 0);
            KRATOS_EXPECT_EQ(return_mapping_variables.Options.Is(FR::PLASTIC_REGION), expected_region != 0);
            KRATOS_EXPECT_EQ(p_flow_rule->GetPlasticRegion(), expected_region);
            KRATOS_EXPECT_EQ(return_mapping_variables.TrialStateFunction > 0.0, expected_region != 0);

            Vector principal_stress = ZeroVector(3);
            for (unsigned int i = 0; i < 3; ++i) {
                principal_stress[i] = stress(i,i);
                for (unsigned int j = 0; j < 3; ++j) {
                    KRATOS_EXPECT_TRUE(std::isfinite(stress(i,j)));
                    KRATOS_EXPECT_NEAR(stress(i,j), i == j ? rExpectedStress[i] : 0.0, 1e-8);
                }
            }
            if (dilatancy_angle == 0.0) {
                KRATOS_EXPECT_NEAR(sum(principal_stress), trial_trace, 1e-8);
            }
            if (expected_region != 0) {
                double yield_function = 0.0;
                p_yield_criterion->CalculateYieldCondition(
                    yield_function, principal_stress, material_properties[COHESION], friction_angle, material_properties);
                KRATOS_EXPECT_NEAR(yield_function, 0.0, 1e-8);
                if (friction_angle == 0.0) {
                    KRATOS_EXPECT_NEAR(principal_stress[0] - principal_stress[2], 2.0 * material_properties[COHESION], 1e-8);
                }
            }

            Matrix tangent = ZeroMatrix(6,6);
            p_flow_rule->ComputeElastoPlasticTangentMatrix(
                return_mapping_variables, strain, 0.0, tangent, material_properties);
            for (unsigned int i = 0; i < 6; ++i) {
                for (unsigned int j = 0; j < 6; ++j) {
                    KRATOS_EXPECT_TRUE(std::isfinite(tangent(i,j)));
                }
            }
            if (friction_angle == 0.0 && expected_region != 0) {
                const double bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio));
                const double shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio));
                const std::array<double,3> yield_normal = {1.0, 0.0, -1.0};
                for (unsigned int i = 0; i < 3; ++i) {
                    for (unsigned int j = 0; j < 3; ++j) {
                        double expected_tangent = bulk_modulus;
                        if (expected_region == 1) {
                            expected_tangent += shear_modulus * (
                                (i == j ? 2.0 : 0.0) - 2.0 / 3.0 - yield_normal[i] * yield_normal[j]);
                        }
                        KRATOS_EXPECT_NEAR(tangent(i,j), expected_tangent, 1e-8);
                    }
                }
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulombZeroFrictionSurface, KratosMPMFastSuite)
    {
        check_mc_return_mapping({20000.0, 5000.0, -10000.0}, {10000.0, 5000.0, 0.0}, 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulombZeroFrictionUpperEdge, KratosMPMFastSuite)
    {
        check_mc_return_mapping({20000.0, 19000.0, -10000.0}, {13000.0, 13000.0, 3000.0}, 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulombZeroFrictionLowerEdge, KratosMPMFastSuite)
    {
        check_mc_return_mapping({20000.0, -9000.0, -10000.0}, {7000.0, -3000.0, -3000.0}, 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulombZeroFrictionElastic, KratosMPMFastSuite)
    {
        check_mc_return_mapping({4000.0, 1000.0, -3000.0}, {4000.0, 1000.0, -3000.0}, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulombApex, KratosMPMFastSuite)
    {
        const double angle = Globals::Pi / 6.0;
        const double apex = 5000.0 / std::tan(angle);
        check_mc_return_mapping({20000.0, 19000.0, 18000.0}, {apex, apex, apex}, 4, angle, angle);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawMohrCoulombSmallFriction, KratosMPMFastSuite)
    {
        for (const double friction_angle : {1e-4, 1e-8, 1e-12, 1e-16}) {
            const double k = (1.0 + std::sin(friction_angle)) / (1.0 - std::sin(friction_angle));
            const double cohesion_coefficient = 10000.0 * std::sqrt(k);

            // Solve the active yield planes together with zero volumetric plastic strain.
            const double delta = (k * 20000.0 + 10000.0 - cohesion_coefficient) / (k + 1.0);
            check_mc_return_mapping(
                {20000.0, 5000.0, -10000.0}, {20000.0 - delta, 5000.0, -10000.0 + delta}, 1, friction_angle);

            const double upper_stress = (29000.0 + cohesion_coefficient) / (k + 2.0);
            check_mc_return_mapping(
                {20000.0, 19000.0, -10000.0}, {upper_stress, upper_stress, k * upper_stress - cohesion_coefficient}, 2, friction_angle);

            const double lower_stress = (1000.0 + 2.0 * cohesion_coefficient) / (1.0 + 2.0 * k);
            check_mc_return_mapping(
                {20000.0, -9000.0, -10000.0}, {lower_stress, k * lower_stress - cohesion_coefficient, k * lower_stress - cohesion_coefficient}, 3, friction_angle);
        }
    }

} // namespace Testing
} // namespace Kratos
