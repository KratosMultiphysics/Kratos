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


namespace Kratos
{
namespace Testing
{

    typedef Node<3> NodeType;

    typedef MCYieldCriterion MC;

    typedef ModifiedCamClayYieldCriterion MCC;

    /**
    * Check the correct yield condition calculation of a given stress condition
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleYieldConditionCalculation, KratosParticleMechanicsFastSuite)
    {

    }

    /**
    * Check the correct calculation of the derivatives of the yield surfaces
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleYieldFunctionFirstDerivativeCalculation, KratosParticleMechanicsFastSuite)
    {

    }

    /**
    * Check the correct calculation of the second derivatives of the yield surfaces
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleYieldFunctionSecondDerivativeCalculation, KratosParticleMechanicsFastSuite)
    {

    }



} // namespace Testing
} // namespace Kratos
