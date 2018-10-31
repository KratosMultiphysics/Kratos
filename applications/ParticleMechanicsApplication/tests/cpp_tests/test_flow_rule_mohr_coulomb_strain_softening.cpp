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
#include "custom_constitutive/hardening_laws/exponential_strain_softening_law.hpp"

// Yields
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"

// Flow Rule
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"


namespace Kratos
{
namespace Testing
{

    typedef Node<3> NodeType;

    typedef MPMHardeningLaw HL;

    typedef ExponentialStrainSofteningLaw ExpSSL;

    typedef MPMYieldCriterion YC;

    typedef MCYieldCriterion MCYC;

    typedef MPMFlowRule FR;

    typedef MCPlasticFlowRule MCFR;

    void GenerateTestMCSSVariables(
        Properties &rMaterialProperties)
    {
        rMaterialProperties.SetValue(YOUNG_MODULUS, 2200.0);
        rMaterialProperties.SetValue(POISSON_RATIO, 2000000.0);
        rMaterialProperties.SetValue(INTERNAL_FRICTION_ANGLE, 0.872664625997165);
        rMaterialProperties.SetValue(INTERNAL_FRICTION_ANGLE_RESIDUAL, 0.5235987755982);
        rMaterialProperties.SetValue(COHESION, 5000.0);
        rMaterialProperties.SetValue(COHESION_RESIDUAL, 0.0);
        rMaterialProperties.SetValue(INTERNAL_DILATANCY_ANGLE, 0.436332312998582);
        rMaterialProperties.SetValue(INTERNAL_DILATANCY_ANGLE_RESIDUAL, 0.0);

    }

    /**
    * Check the correct calculation of the integrated stress with the CL's
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleConstitutiveLawMohrCoulombStrainSoftening, KratosParticleMechanicsFastSuite)
    {
        Properties material_properties;
        GenerateTestMCSSVariables(material_properties);

        // Construct Flow Rule
        HL::Pointer mcss_hl_pointer = HL::Pointer( new ExpSSL() );
        YC::Pointer mcss_yc_pointer = YC::Pointer( new MCYC() );
        FR::Pointer mcss_fr_pointer = FR::Pointer( new MCFR() );
        mcss_fr_pointer->InitializeMaterial(mcss_yc_pointer, mcss_hl_pointer, material_properties);



    }

} // namespace Testing
} // namespace Kratos
