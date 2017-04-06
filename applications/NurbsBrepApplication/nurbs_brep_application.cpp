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
  
  KRATOS_REGISTER_VARIABLE( LOCAL_PARAMETERS)
  KRATOS_REGISTER_VARIABLE( FACE_BREP_ID)

  KRATOS_REGISTER_VARIABLE( CONTROL_POINT_IDS)

  KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_VALUES)
  KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_DERIVATIVES)
  KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_SECOND_DERIVATIVES)
  //KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
  //KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
  
  //// edge integration
  //KRATOS_REGISTER_VARIABLE( TANGENTS)
  //
  //// penalty factor
  //KRATOS_REGISTER_VARIABLE( PENALTY_FACTOR)

  //// coupling and support
  //KRATOS_REGISTER_VARIABLE(DISPLACEMENT_ROTATION_FIX)
  //// for load condition
  //KRATOS_REGISTER_VARIABLE(LOAD_TYPE)
  //KRATOS_REGISTER_VARIABLE( DISTRIBUTED_LOAD_FACTOR)
}
}  // namespace Kratos.
