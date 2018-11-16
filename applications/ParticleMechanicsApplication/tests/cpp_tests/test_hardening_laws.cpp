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

// Hardening laws
#include "custom_constitutive/hardening_laws/cam_clay_hardening_law.hpp"
#include "custom_constitutive/hardening_laws/exponential_strain_softening_law.hpp"


namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    typedef Node<3> NodeType;

    typedef MPMHardeningLaw HL;

    typedef CamClayHardeningLaw CCHL;

    typedef ExponentialStrainSofteningLaw ExpSSHL;

    void GenerateTestVariables(
        Properties &rMaterialProperties)
    {
        rMaterialProperties.SetValue(SWELLING_SLOPE, 0.018);
        rMaterialProperties.SetValue(NORMAL_COMPRESSION_SLOPE, 0.130);

        rMaterialProperties.SetValue(COHESION, 5000.0);
        rMaterialProperties.SetValue(COHESION_RESIDUAL, 0.0);
        rMaterialProperties.SetValue(INTERNAL_FRICTION_ANGLE, 0.872664625997165);
        rMaterialProperties.SetValue(INTERNAL_FRICTION_ANGLE_RESIDUAL, 0.5235987755982);
        rMaterialProperties.SetValue(INTERNAL_DILATANCY_ANGLE, 0.436332312998582);
        rMaterialProperties.SetValue(INTERNAL_DILATANCY_ANGLE_RESIDUAL, 0.0);
        rMaterialProperties.SetValue(SHAPE_FUNCTION_BETA, 4.0);

    }

    /**
    * Check the correct hardening parameter calculation of a given stress condition
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleHardeningParameterCalculation, KratosParticleMechanicsFastSuite)
    {
        Properties material_properties;
        GenerateTestVariables(material_properties);

        // Analytical solutions of the yield surface
        double CC_analytical;
        CC_analytical  = -48892.40454;
        Vector ExpSS_analytical = ZeroVector(3);
        ExpSS_analytical[0] = -0.9359435085; // Internal friction angle
        ExpSS_analytical[1] = -13406.40092; // Cohesion
        ExpSS_analytical[2] = -1.169929185; // Internal dilation angle

        // Construct Hardening law
        // 1. For Cam Clay hardening law
        HL::Pointer cc_hl_pointer = HL::Pointer( new CCHL() );
        cc_hl_pointer->InitializeMaterial(material_properties);

        // 2. For Exponential strain softening law
        HL::Pointer ExpSS_hl_pointer = HL::Pointer( new ExpSSHL() );
        ExpSS_hl_pointer->InitializeMaterial(material_properties);

        // Test parameters
        double CC_result = 0.0;
        Vector ExpSS_result = ZeroVector(3);

        double alpha = 0.1;
        double prev_pc = -119.4e3;

        CC_result = cc_hl_pointer->CalculateHardening(CC_result, alpha, prev_pc);
        ExpSS_result[0] = ExpSS_hl_pointer->CalculateHardening(ExpSS_result[0], alpha, INTERNAL_FRICTION_ANGLE);
        ExpSS_result[1] = ExpSS_hl_pointer->CalculateHardening(ExpSS_result[1], alpha, COHESION);
        ExpSS_result[2] = ExpSS_hl_pointer->CalculateHardening(ExpSS_result[2], alpha, INTERNAL_DILATANCY_ANGLE);

        // Check results
        KRATOS_CHECK_NEAR(CC_result, CC_analytical, 1.e-3);
        KRATOS_CHECK_NEAR(ExpSS_result[0], ExpSS_analytical[0], tolerance);
        KRATOS_CHECK_NEAR(ExpSS_result[1], ExpSS_analytical[1], tolerance);
        KRATOS_CHECK_NEAR(ExpSS_result[2], ExpSS_analytical[2], tolerance);

    }

} // namespace Testing
} // namespace Kratos
