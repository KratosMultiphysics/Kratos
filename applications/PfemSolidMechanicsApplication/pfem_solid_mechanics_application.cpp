//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
  //Create Variables

  //solution
  KRATOS_CREATE_VARIABLE(int, NUMBER_OF_ACTIVE_CONTACTS );
  KRATOS_CREATE_VARIABLE(int, NUMBER_OF_STICK_CONTACTS );
  KRATOS_CREATE_VARIABLE(int, NUMBER_OF_SLIP_CONTACTS );
  
  //constitutive law	
  KRATOS_CREATE_VARIABLE(double, MEAN_ERROR );

  //element

  //thermal

  //mechanical
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_CONTACT_NORMAL );
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_CONTACT_TANGENT );

  //geometrical
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( OFFSET );
  KRATOS_CREATE_VARIABLE(Vector, BOUNDARY_NORMAL );
  KRATOS_CREATE_VARIABLE(double, SHRINK_FACTOR );

  //domain definition
  KRATOS_CREATE_VARIABLE(unsigned int, DOMAIN_LABEL );
  KRATOS_CREATE_VARIABLE(bool        , RIGID_WALL);
  KRATOS_CREATE_VARIABLE(double      , WALL_TIP_RADIUS );
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT );
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY );

  //contact condition
  KRATOS_CREATE_VARIABLE(Condition::Pointer, MASTER_CONDITION );
  KRATOS_CREATE_VARIABLE(WeakPointerVector< Element >, MASTER_ELEMENTS );
  KRATOS_CREATE_VARIABLE(WeakPointerVector< Node<3> >, MASTER_NODES );

  //properties
  KRATOS_CREATE_VARIABLE(bool, PENALTY_CONTACT );
  KRATOS_CREATE_VARIABLE(bool, FRICTION_ACTIVE );
  KRATOS_CREATE_VARIABLE(double, TAU_STAB );
  KRATOS_CREATE_VARIABLE(double, MU_STATIC );
  KRATOS_CREATE_VARIABLE(double, MU_DYNAMIC );

  //global flags
  // KRATOS_CREATE_FLAG( FLUID,         16 );
  // KRATOS_CREATE_FLAG( STRUCTURE,     17 );
  // KRATOS_CREATE_FLAG( SOLID,         18 );
  // KRATOS_CREATE_FLAG( RIGID,         19 );
  // KRATOS_CREATE_FLAG( CONTACT,       20 );
  
  
  // KRATOS_CREATE_FLAG( BOUNDARY,      21 );
  // KRATOS_CREATE_FLAG( FREE_SURFACE,  22 );
  
  // KRATOS_CREATE_FLAG( INTERFACE,     23 );
  
  // KRATOS_CREATE_FLAG( ENGAGED,       24 );
  // KRATOS_CREATE_FLAG( ISOLATED,      25 );
  
  // KRATOS_CREATE_FLAG( REFINE,        26 );
  // KRATOS_CREATE_FLAG( INSERTED,      27 );
  // KRATOS_CREATE_FLAG( RELEASE,       28 );

  KratosPfemSolidMechanicsApplication::KratosPfemSolidMechanicsApplication():
    mCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSkinMultipleCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mSkinMultipleCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mWallTipCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mWallTipCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mContactDomain2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymContactDomain2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    
  {}
  
  void KratosPfemSolidMechanicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosPfemSolidMechanicsApplication... " << std::endl;
    
    //Register Elements

    //Register Conditions

    KRATOS_REGISTER_CONDITION( "Condition2D", mCondition2D );
    KRATOS_REGISTER_CONDITION( "Condition3D", mCondition3D );

    KRATOS_REGISTER_CONDITION( "SkinCondition2D", mSkinMultipleCondition2D );
    KRATOS_REGISTER_CONDITION( "SkinCondition3D", mSkinMultipleCondition3D );
   
    KRATOS_REGISTER_CONDITION( "WallTipCondition2D", mWallTipCondition2D );
    KRATOS_REGISTER_CONDITION( "WallTipCondition3D", mWallTipCondition3D );

    KRATOS_REGISTER_CONDITION( "ContactDomain2DCondition", mContactDomain2DCondition );

    KRATOS_REGISTER_CONDITION( "AxisymContactDomain2DCondition",   mAxisymContactDomain2DCondition );

    //Register Constitutive Laws
    
    //Register Variables

    //solution
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_ACTIVE_CONTACTS );
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_STICK_CONTACTS );
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_SLIP_CONTACTS );
    
    //constitutive law	
    KRATOS_REGISTER_VARIABLE( MEAN_ERROR );

    //element

    //thermal

    //mechanical
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FORCE_CONTACT_NORMAL );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FORCE_CONTACT_TANGENT );

    //geometrical
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( OFFSET );
    KRATOS_REGISTER_VARIABLE( BOUNDARY_NORMAL );
    KRATOS_REGISTER_VARIABLE( SHRINK_FACTOR );

    //domain definition
    KRATOS_REGISTER_VARIABLE( DOMAIN_LABEL );
    KRATOS_REGISTER_VARIABLE( RIGID_WALL );
    KRATOS_REGISTER_VARIABLE( WALL_TIP_RADIUS );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY );

    //contact condition
    KRATOS_REGISTER_VARIABLE( MASTER_CONDITION );
    KRATOS_REGISTER_VARIABLE( MASTER_ELEMENTS );
    KRATOS_REGISTER_VARIABLE( MASTER_NODES );

    //properties
    KRATOS_REGISTER_VARIABLE( PENALTY_CONTACT );
    KRATOS_REGISTER_VARIABLE( FRICTION_ACTIVE );
    KRATOS_REGISTER_VARIABLE( TAU_STAB );
    KRATOS_REGISTER_VARIABLE( MU_STATIC );
    KRATOS_REGISTER_VARIABLE( MU_DYNAMIC );
  
    //flags
    // KRATOS_REGISTER_FLAG( FLUID );
    // KRATOS_REGISTER_FLAG( STRUCTURE );
    // KRATOS_REGISTER_FLAG( SOLID );
    // KRATOS_REGISTER_FLAG( RIGID );
    // KRATOS_REGISTER_FLAG( CONTACT );

    // KRATOS_REGISTER_FLAG( BOUNDARY );
    // KRATOS_REGISTER_FLAG( FREE_SURFACE );

    // KRATOS_REGISTER_FLAG( INTERFACE );

    // KRATOS_REGISTER_FLAG( ENGAGED );
    // KRATOS_REGISTER_FLAG( ISOLATED );

    // KRATOS_REGISTER_FLAG( REFINE );
    // KRATOS_REGISTER_FLAG( INSERTED );
    // KRATOS_REGISTER_FLAG( RELEASE );      
  }
  
}  // namespace Kratos.


