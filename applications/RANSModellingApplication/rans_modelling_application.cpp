//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//


// System includes


// External includes


// Project includes
#include "rans_modelling_application.h"
#include "rans_modelling_application_variables.h"


namespace Kratos {

KratosRANSModellingApplication::KratosRANSModellingApplication():
    KratosApplication("RANSModellingApplication")
    {}

void KratosRANSModellingApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosRANSModellingApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY )
  KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE )
  KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_RATE )
  KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_2 )
  KRATOS_REGISTER_VARIABLE( IS_CO_SOLVING_PROCESS_ACTIVE )
  KRATOS_REGISTER_VARIABLE( OLD_CONVERGENCE_VARIABLE )
  KRATOS_REGISTER_VARIABLE( RANS_Y_PLUS )
  KRATOS_REGISTER_VARIABLE( RANS_AUXILIARY_VARIABLE_1 )
  KRATOS_REGISTER_VARIABLE( RANS_AUXILIARY_VARIABLE_2 )
  KRATOS_REGISTER_VARIABLE( WALL_SMOOTHNESS_BETA )
  KRATOS_REGISTER_VARIABLE( WALL_VON_KARMAN )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C_MU )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C1 )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C2 )
  KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY_MIN )
  KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY_MAX )
  KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA )
  KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )

}
}  // namespace Kratos.
