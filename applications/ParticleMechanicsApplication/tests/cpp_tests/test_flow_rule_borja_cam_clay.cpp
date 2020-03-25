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

// External includes

// Project includes
#include "testing/testing.h"

// Application includes
#include "particle_mechanics_application_variables.h"

// Hardening law
#include "custom_constitutive/hardening_laws/cam_clay_hardening_law.hpp"

// Yields
#include "custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.hpp"

// Flow Rule
#include "custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.hpp"


namespace Kratos
{
namespace Testing
{
    typedef Node<3> NodeType;

    typedef ParticleHardeningLaw HL;

    typedef CamClayHardeningLaw CCHL;

    typedef ParticleYieldCriterion YC;

    typedef ModifiedCamClayYieldCriterion MCCYC;

    typedef ParticleFlowRule FR;

    typedef BorjaCamClayPlasticFlowRule BCCFR;

    void GenerateTestBCCVariables(
        Matrix& rStress, Matrix& rStrain,
        Properties &rMaterialProperties)
    {
        rStress = ZeroMatrix(3,3);
        rStrain = ZeroMatrix(3,3);
        rStrain(0,0) = 5.5e-2;
        rStrain(1,1) = -12.2e-2;
        rStrain(2,2) = 8.3e-2;

        rMaterialProperties.SetValue(PRE_CONSOLIDATION_STRESS, -90000.0);
        rMaterialProperties.SetValue(OVER_CONSOLIDATION_RATIO, 1.0);
        rMaterialProperties.SetValue(SWELLING_SLOPE          , 0.018);
        rMaterialProperties.SetValue(NORMAL_COMPRESSION_SLOPE, 0.13);
        rMaterialProperties.SetValue(CRITICAL_STATE_LINE     , 1.05);
        rMaterialProperties.SetValue(INITIAL_SHEAR_MODULUS   , 5400000.0);
        rMaterialProperties.SetValue(ALPHA_SHEAR             , 0.0);

    }

    /**
    * Check the correct calculation of the integrated stress with the CL's
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleConstitutiveLawBorjaCamClay, KratosParticleMechanicsFastSuite)
    {
        Matrix stress, strain;
        Properties material_properties;
        GenerateTestBCCVariables(stress, strain, material_properties);

        // Construct Flow Rule
        HL::Pointer bcc_hl_pointer = HL::Pointer( new CCHL() );
        YC::Pointer bcc_yc_pointer = YC::Pointer( new MCCYC() );
        FR::Pointer bcc_fr_pointer = FR::Pointer( new BCCFR() );
        bcc_fr_pointer->InitializeMaterial(bcc_yc_pointer, bcc_hl_pointer, material_properties);

        // Compute trial elastic stresses
        FR::RadialReturnVariables rma_variables;
        bcc_fr_pointer->CalculatePrincipalStressTrial(rma_variables, strain, stress);

        Vector stress_trial_analytic = ZeroVector(3);
        stress_trial_analytic[0] =  4.993998938544e+05;
        stress_trial_analytic[1] = -1.412200106146e+06;
        stress_trial_analytic[2] =  8.017998938544e+05;

        KRATOS_CHECK_NEAR(stress(0,0), stress_trial_analytic[0], 1e-6);
        KRATOS_CHECK_NEAR(stress(1,1), stress_trial_analytic[1], 1e-6);
        KRATOS_CHECK_NEAR(stress(2,2), stress_trial_analytic[2], 1e-6);

        // Compute new stresses after return mapping
        Matrix dummy_deformation_gradient = IdentityMatrix(3);
        bcc_fr_pointer->CalculateReturnMapping( rma_variables, dummy_deformation_gradient, stress, strain);

        Vector stress_analytic = ZeroVector(3);
        stress_analytic[0] = -3.112443504093e+04;
        stress_analytic[1] = -7.351387390960e+04;
        stress_analytic[2] = -2.441876109561e+04;
        const double yield_analytic = 3.919787896261e+12;

        KRATOS_CHECK_NEAR(rma_variables.TrialStateFunction, yield_analytic, 1.0);
        KRATOS_CHECK_NEAR(stress(0,0), stress_analytic[0], 1e-6);
        KRATOS_CHECK_NEAR(stress(1,1), stress_analytic[1], 1e-6);
        KRATOS_CHECK_NEAR(stress(2,2), stress_analytic[2], 1e-6);

    }

} // namespace Testing
} // namespace Kratos
