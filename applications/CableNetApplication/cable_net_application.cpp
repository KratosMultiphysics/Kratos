//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "cable_net_application.h"
#include "cable_net_application_variables.h"


namespace Kratos {

KratosCableNetApplication::KratosCableNetApplication():
    KratosApplication("CableNetApplication")
    {}

void KratosCableNetApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosCableNetApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( DOF_1 )
  KRATOS_REGISTER_VARIABLE( DOF_2 )
  KRATOS_REGISTER_VARIABLE( ScalarVariable )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VectorVariable )

}
}  // namespace Kratos.
