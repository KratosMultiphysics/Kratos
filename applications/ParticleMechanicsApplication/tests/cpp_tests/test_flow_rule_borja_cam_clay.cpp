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

    typedef MPMHardeningLaw HL;

    typedef CamClayHardeningLaw CCHL;

    typedef MPMYieldCriterion YC;

    typedef ModifiedCamClayYieldCriterion MCCYC;

    typedef MPMFlowRule FR;

    typedef BorjaCamClayPlasticFlowRule BCCFR;

    void GenerateTestBCCVariables(
        Properties &rMaterialProperties)
    {
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
        Properties material_properties;
        GenerateTestBCCVariables(material_properties);

        // Construct Flow Rule
        HL::Pointer bcc_hl_pointer = HL::Pointer( new CCHL() );
        YC::Pointer bcc_yc_pointer = YC::Pointer( new MCCYC() );
        FR::Pointer bcc_fr_pointer = FR::Pointer( new BCCFR() );
        bcc_fr_pointer->InitializeMaterial(bcc_yc_pointer, bcc_hl_pointer, material_properties);



    }

} // namespace Testing
} // namespace Kratos
