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

    typedef Node<3> NodeType;

    typedef CamClayHardeningLaw CCH;

    typedef ExponentialStrainSofteningLaw ExpSS;

    /**
    * Check the correct hardening parameter calculation of a given stress condition
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleHardeningParameterCalculation, KratosParticleMechanicsFastSuite)
    {

    }

} // namespace Testing
} // namespace Kratos
