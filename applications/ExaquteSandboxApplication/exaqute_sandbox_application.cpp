//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes


// External includes


// Project includes
#include "exaqute_sandbox_application.h"
#include "exaqute_sandbox_application_variables.h"


namespace Kratos {

KratosExaquteSandboxApplication::KratosExaquteSandboxApplication():
    KratosApplication("ExaquteSandboxApplication")
    {}

void KratosExaquteSandboxApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosExaquteSandboxApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( AVERAGED_DIVERGENCE )
  KRATOS_REGISTER_VARIABLE( DIVERGENCE )
  KRATOS_REGISTER_VARIABLE( VELOCITY_H1_SEMINORM )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_NOISE )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AVERAGED_VELOCITY )
  KRATOS_REGISTER_VARIABLE( AVERAGED_PRESSURE )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DRAG_FORCE_CENTER_AUX)

}
}  // namespace Kratos.
