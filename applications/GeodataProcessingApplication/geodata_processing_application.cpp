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
#include "geodata_processing_application.h"
#include "geodata_processing_application_variables.h"


namespace Kratos {

KratosGeodataProcessingApplication::KratosGeodataProcessingApplication():
    KratosApplication("GeodataProcessingApplication")
    {}

void KratosGeodataProcessingApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosGeodataProcessingApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( EXTRUSION_HEIGHT )
  KRATOS_REGISTER_VARIABLE( DOF_2 )
  KRATOS_REGISTER_VARIABLE( ScalarVariable )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VectorVariable )

}
}  // namespace Kratos.
