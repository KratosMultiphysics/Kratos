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

// Yields
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
#include "custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.hpp"

// Hardening law
#include "custom_constitutive/hardening_laws/MPM_hardening_law.hpp"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    const double GetPI()
    {
        return std::atan(1.0)*4.0;
    }

    typedef MPMHardeningLaw HL;

    typedef CamClayHardeningLaw CCHL;

    typedef MPMYieldCriterion YC;

    typedef MCYieldCriterion MC;

    typedef ModifiedCamClayYieldCriterion MCC;

    void GenerateTestVariables(
        Vector &rStressVector,
        Properties &rMaterialProperties)
    {
        rStressVector = ZeroVector(3);
        rStressVector[0] = -5.0e6;
        rStressVector[1] = -6.0e6;
        rStressVector[2] = -1.0e7;

        rMaterialProperties.SetValue(CRITICAL_STATE_LINE, 1.05);
        rMaterialProperties.SetValue(SWELLING_SLOPE, 0.018);
        rMaterialProperties.SetValue(NORMAL_COMPRESSION_SLOPE, 0.130);

    }

    /**
    * Check the correct yield condition calculation of a given stress condition
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleYieldConditionCalculation, KratosParticleMechanicsFastSuite)
    {
        Vector principal_stress = ZeroVector(3);
        Properties material_properties;
        GenerateTestVariables(principal_stress, material_properties);

        // Analytical solutions of the yield surface
        double MC_analytical, MCC_analytical;
        MC_analytical = -5017320.508;
        MCC_analytical = 6.770537222e13;

        // Construct YieldCriterion
        // 1. For Mohr Coulomb yield criterion
        HL::Pointer mc_hl_pointer = HL::Pointer( new HL() );
        YC::Pointer mc_yc_pointer = YC::Pointer( new MC());
        mc_yc_pointer->InitializeMaterial(mc_hl_pointer, material_properties);

        // 2. For Modified Cam Clay yield criterion
        HL::Pointer mcc_hl_pointer = HL::Pointer( new CCHL() );
        YC::Pointer mcc_yc_pointer = YC::Pointer( new MCC());
        mcc_yc_pointer->InitializeMaterial(mcc_hl_pointer, material_properties);

        // Test parameters
        double MC_result = 0.0, MCC_result = 0.0;
        double cohesion = 5000.0;
        double friction_angle = GetPI()/6;
        double alpha = 0.1;
        double prev_pc = -119.4e3;

        MC_result  = mc_yc_pointer->CalculateYieldCondition(MC_result, principal_stress, cohesion, friction_angle);
        MCC_result = mcc_yc_pointer->CalculateYieldCondition(MCC_result, principal_stress, alpha, prev_pc);

        // Check results
        KRATOS_CHECK_NEAR(MC_result, MC_analytical, 1);
        KRATOS_CHECK_NEAR(MCC_result, MCC_analytical, 1.e5);

    }

    /**
    * Check the correct calculation of the derivatives of the yield surfaces
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleYieldFunctionFirstDerivativeCalculation, KratosParticleMechanicsFastSuite)
    {
        Vector principal_stress = ZeroVector(3);
        Properties material_properties;
        GenerateTestVariables(principal_stress, material_properties);

        // Analytical solutions of the yield surface
        Vector MCC_derivative_analytical = ZeroVector(3);
        MCC_derivative_analytical[0] = -13951107.60;
        MCC_derivative_analytical[1] = 8313062.485;
        MCC_derivative_analytical[2] = 7000000.0;

        // Construct YieldCriterion
        // 1. For Modified Cam Clay yield criterion
        HL::Pointer mcc_hl_pointer = HL::Pointer( new CCHL() );
        YC::Pointer mcc_yc_pointer = YC::Pointer( new MCC());
        mcc_yc_pointer->InitializeMaterial(mcc_hl_pointer, material_properties);

        // Test parameters
        Vector MCC_derivative_result = ZeroVector(3);
        double alpha = 0.1;
        double prev_pc = -119.4e3;

        mcc_yc_pointer->CalculateYieldFunctionDerivative(principal_stress, MCC_derivative_result, alpha, prev_pc);

        // Check results
        KRATOS_CHECK_NEAR(MCC_derivative_result[0], MCC_derivative_analytical[0], 1);
        KRATOS_CHECK_NEAR(MCC_derivative_result[1], MCC_derivative_analytical[1], 1);
        KRATOS_CHECK_NEAR(MCC_derivative_result[2], MCC_derivative_analytical[2], 1);

    }

    /**
    * Check the correct calculation of the second derivatives of the yield surfaces
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleYieldFunctionSecondDerivativeCalculation, KratosParticleMechanicsFastSuite)
    {
        Vector principal_stress = ZeroVector(3);
        Properties material_properties;
        GenerateTestVariables(principal_stress, material_properties);

        // Analytical solutions of the yield surface
        Vector MCC_derivative2_analytical = ZeroVector(6);
        MCC_derivative2_analytical[0] = 2.0;
        MCC_derivative2_analytical[1] = 1.814059;
        MCC_derivative2_analytical[2] = 0.0;
        MCC_derivative2_analytical[3] = 0.0;
        MCC_derivative2_analytical[4] = 0.0;
        MCC_derivative2_analytical[5] =-1.0;

        // Construct YieldCriterion
        // 1. For Modified Cam Clay yield criterion
        HL::Pointer mcc_hl_pointer = HL::Pointer( new CCHL() );
        YC::Pointer mcc_yc_pointer = YC::Pointer( new MCC());
        mcc_yc_pointer->InitializeMaterial(mcc_hl_pointer, material_properties);

        // Test parameters
        Vector MCC_derivative2_result = ZeroVector(6);

        mcc_yc_pointer->CalculateYieldFunctionSecondDerivative(principal_stress, MCC_derivative2_result);

        // Check results
        KRATOS_CHECK_NEAR(MCC_derivative2_result[0], MCC_derivative2_analytical[0], tolerance);
        KRATOS_CHECK_NEAR(MCC_derivative2_result[1], MCC_derivative2_analytical[1], tolerance);
        KRATOS_CHECK_NEAR(MCC_derivative2_result[2], MCC_derivative2_analytical[2], tolerance);
        KRATOS_CHECK_NEAR(MCC_derivative2_result[3], MCC_derivative2_analytical[3], tolerance);
        KRATOS_CHECK_NEAR(MCC_derivative2_result[4], MCC_derivative2_analytical[4], tolerance);
        KRATOS_CHECK_NEAR(MCC_derivative2_result[5], MCC_derivative2_analytical[5], tolerance);

    }




} // namespace Testing
} // namespace Kratos
