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
#include "includes/process_info.h"
#include "testing/testing.h"

// Application includes
#include "particle_mechanics_application_variables.h"

// Hardening law
#include "custom_constitutive/hardening_laws/MPM_hardening_law.hpp"

// Yields
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"

// Flow Rule
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"


namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    typedef Node<3> NodeType;

    typedef MPMHardeningLaw HL;

    typedef MPMYieldCriterion YC;

    typedef MCYieldCriterion MCYC;

    typedef MPMFlowRule FR;

    typedef MCPlasticFlowRule MCFR;

    void GenerateTestMCVariables(
        Matrix& rStress, Matrix& rStrain,
        Properties &rMaterialProperties)
    {
        rStress = ZeroMatrix(3);
        rStrain = ZeroMatrix(3);
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
    KRATOS_TEST_CASE_IN_SUITE(ParticleConstitutiveLawMohrCoulomb, KratosParticleMechanicsFastSuite)
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
        mc_fr_pointer->CalculatePrincipalStressTrial(rma_variables, strain, stress);

        Vector stress_trial_analytic = ZeroVector(3);
        stress_trial_analytic[0] =  1.03076923077e+05;
        stress_trial_analytic[1] = -1.69230769230e+05;
        stress_trial_analytic[2] =  1.46153846154e+05;

        KRATOS_CHECK_NEAR(stress(0,0), stress_trial_analytic[0], tolerance);
        KRATOS_CHECK_NEAR(stress(1,1), stress_trial_analytic[1], tolerance);
        KRATOS_CHECK_NEAR(stress(2,2), stress_trial_analytic[2], tolerance);

        // Compute new stresses after return mapping
        Matrix dummy_deformation_gradient = IdentityMatrix(3);
        mc_fr_pointer->CalculateReturnMapping( rma_variables, dummy_deformation_gradient, stress, strain);

        Vector stress_analytic = ZeroVector(3);
        stress_analytic[0] = -5.31271630223e+04;
        stress_analytic[1] = -1.76701997143e+05;
        stress_analytic[2] = -5.31271630223e+04;
        const double yield_analytic = 5.90371799617e+05;

        KRATOS_CHECK_NEAR(rma_variables.TrialStateFunction, yield_analytic, tolerance);
        KRATOS_CHECK_NEAR(stress(0,0), stress_analytic[0], tolerance);
        KRATOS_CHECK_NEAR(stress(1,1), stress_analytic[1], tolerance);
        KRATOS_CHECK_NEAR(stress(2,2), stress_analytic[2], tolerance);

    }

} // namespace Testing
} // namespace Kratos
