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

// Constitutive law
#include "custom_constitutive/hencky_mc_strain_softening_3D_law.hpp"


namespace Kratos
{
namespace Testing
{

    typedef Node<3> NodeType;

    /**
    * Check the correct calculation of the integrated stress with the CL's
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleConstitutiveLawMohrCoulombStrainSoftening, KratosParticleMechanicsFastSuite)
    {

    }

} // namespace Testing
} // namespace Kratos
