//
//   Project Name:        KratosPfemApplication     $
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
#include "pfem_application.h"

namespace Kratos
{
  //Create Variables


  KratosPfemApplication    ::KratosPfemApplication    ():
    mCompositeCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mCompositeCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) )    
  {}
  
  void KratosPfemApplication    ::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
       
    std::cout << "            ___  __                     " << std::endl;
    std::cout << "     KRATOS| _ \\/ _|___ _ __            " << std::endl;
    std::cout << "           |  _/  _/ -_) '  \\            " << std::endl;
    std::cout << "           |_| |_| \\___|_|_|_|APPLICATION " << std::endl;
    std::cout << "Initializing KratosPfemApplication    ...        " << std::endl;                                

    //Register Variables (variables created in pfem_solid_mechanics_application_variables.cpp)
    
    //geometrical definition
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
    KRATOS_REGISTER_VARIABLE( SHRINK_FACTOR )


    //domain definition
    KRATOS_REGISTER_VARIABLE( INITIALIZED_DOMAINS )
    KRATOS_REGISTER_VARIABLE( MESHING_STEP_PERFORMED )
    KRATOS_REGISTER_VARIABLE( MODEL_PART_NAME )

    //boundary definition
    KRATOS_REGISTER_VARIABLE( RIGID_WALL )
    KRATOS_REGISTER_VARIABLE( MASTER_CONDITION )
    KRATOS_REGISTER_VARIABLE( MASTER_ELEMENTS )
    KRATOS_REGISTER_VARIABLE( MASTER_NODES )

    //condition variables
    KRATOS_REGISTER_VARIABLE( CHILDREN_CONDITIONS )

    //modeler criteria
    KRATOS_REGISTER_VARIABLE( MEAN_ERROR )

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "CompositeCondition2D2N", mCompositeCondition2D2N )
    KRATOS_REGISTER_CONDITION( "CompositeCondition3D3N", mCompositeCondition3D3N )



  }
  
}  // namespace Kratos.


