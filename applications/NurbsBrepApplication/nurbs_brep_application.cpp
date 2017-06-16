//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/NURBSBRepApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"  

#include "geometries/geometry.h"

//#include "custom_geometries/meshless_geometry.h"


namespace Kratos {

KratosNurbsBrepApplication::KratosNurbsBrepApplication(){}

void KratosNurbsBrepApplication::Register() {
   // calling base class register to register Kratos components
   KratosApplication::Register();
   //std::cout << "Kratos" << std::endl;
   std::cout << "    _______   ____ _______________________  _________ _________________________________________" << std::endl;
   std::cout << "    \\      \\ |    |   \\______   \\______   \\/  _____ / \\______   \\______   \\_   _____/ \\______   \\" << std::endl;
   std::cout << "    /   |   \\|    |   /|       _/|     | _/\\_____  \\   |    |  _/|       _/|    __)_  |     ___ /" << std::endl;
   std::cout << "   /    |    \\    |  / |    |   \\|     |  \\/        \\  |    |   \\|    |   \\|        \\ |    |" << std::endl;
   std::cout << "   \\____|__  /______/  |____|_  /|______  /_______  /  |______  /|____|_  /_______  / |____|" << std::endl;
   std::cout << "           \\/                 \\/        \\/        \\/          \\/        \\/        \\/" << std::endl;
   std::cout << "Initializing KratosNurbsBrepApplication... " << std::endl;

  //
  KRATOS_REGISTER_VARIABLE( CONTROL_POINT_WEIGHT)

  KRATOS_REGISTER_VARIABLE( INTEGRATION_WEIGHT)

  KRATOS_REGISTER_VARIABLE( TANGENTS_BASIS_VECTOR)
  KRATOS_REGISTER_VARIABLE( LOCAL_PARAMETERS)
  KRATOS_REGISTER_VARIABLE( FACE_BREP_ID)

  KRATOS_REGISTER_VARIABLE( CONTROL_POINT_IDS)

  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_VALUES)
  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_DERIVATIVES)
  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_SECOND_DERIVATIVES)

  //Slave Variables
  KRATOS_REGISTER_VARIABLE( LOCATION_SLAVE)
  KRATOS_REGISTER_VARIABLE( LOCAL_PARAMETERS_SLAVE)

  KRATOS_REGISTER_VARIABLE( TANGENTS_BASIS_VECTOR_SLAVE)
  KRATOS_REGISTER_VARIABLE( FACE_BREP_ID_SLAVE)

  KRATOS_REGISTER_VARIABLE( CONTROL_POINT_IDS_SLAVE)

  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_SLAVE)
  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_DERIVATIVES_SLAVE)
  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_SECOND_DERIVATIVES_SLAVE)
}
}  // namespace Kratos.
