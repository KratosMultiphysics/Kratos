//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes 

// Project includes
#include "includes/define.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/line_2d.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

// Core applications
#include "pfem_base_application.h"

namespace Kratos
{
  //Create Variables


  KratosPfemBaseApplication::KratosPfemBaseApplication()    
  {}
  
  void KratosPfemBaseApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
       
    std::cout << "            ___  __           ___                           " << std::endl;
    std::cout << "     KRATOS| _ \\/ _|___ _ __ | _ ) __ _ ___ ___             " << std::endl;
    std::cout << "           |  _/  _/ -_) '  \\| _ \\/ _` (_-</ -_)            " << std::endl;
    std::cout << "           |_| |_| \\___|_|_|_|___/\\__,_/__/\\___|APPLICATION " << std::endl;
    std::cout << "Initializing KratosPfemBaseApplication...                   " << std::endl;                                

    //Register Variables (variables created in pfem_solid_mechanics_application_variables.cpp)
    
    //geometrical definition
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
    KRATOS_REGISTER_VARIABLE( SHRINK_FACTOR )
    //KRATOS_REGISTER_VARIABLE( NORMAL )
    //KRATOS_REGISTER_VARIABLE( NODAL_H )

    //domain definition
    KRATOS_REGISTER_VARIABLE( DOMAIN_LABEL )
    //KRATOS_REGISTER_VARIABLE( RIGID_WALL )
    //KRATOS_REGISTER_VARIABLE( WALL_TIP_RADIUS )
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT )
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY )

    //modeler criteria
    KRATOS_REGISTER_VARIABLE( MEAN_ERROR )

    //Register Conditions


  }
  
}  // namespace Kratos.


