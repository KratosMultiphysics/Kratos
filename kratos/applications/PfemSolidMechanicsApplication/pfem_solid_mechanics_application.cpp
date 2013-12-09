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
  KRATOS_CREATE_VARIABLE(int, NUMBER_OF_ACTIVE_CONTACTS )
  KRATOS_CREATE_VARIABLE(int, NUMBER_OF_STICK_CONTACTS )
  KRATOS_CREATE_VARIABLE(int, NUMBER_OF_SLIP_CONTACTS )
  
  //constitutive law	
  KRATOS_CREATE_VARIABLE(double, MEAN_ERROR )

  //element

  //thermal

  //mechanical

  //geometrical
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
  KRATOS_CREATE_VARIABLE(Vector, BOUNDARY_NORMAL )
  KRATOS_CREATE_VARIABLE(double, SHRINK_FACTOR )

  //domain definition
  KRATOS_CREATE_VARIABLE(unsigned int, DOMAIN_LABEL )
  KRATOS_CREATE_VARIABLE(int         , RIGID_WALL )
  KRATOS_CREATE_VARIABLE(double      , WALL_TIP_RADIUS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY )

  //contact condition
  KRATOS_CREATE_VARIABLE(Condition::Pointer, MASTER_CONDITION )
  KRATOS_CREATE_VARIABLE(WeakPointerVector< Element >, MASTER_ELEMENTS )
  KRATOS_CREATE_VARIABLE(WeakPointerVector< Node<3> >, MASTER_NODES )

  //contact 
  KRATOS_CREATE_VARIABLE(bool, FRICTION_ACTIVE )
  KRATOS_CREATE_VARIABLE(double, PENALTY_PARAMETER )
  KRATOS_CREATE_VARIABLE(double, TAU_STAB )
  KRATOS_CREATE_VARIABLE(double, MU_STATIC )
  KRATOS_CREATE_VARIABLE(double, MU_DYNAMIC )

  KratosPfemSolidMechanicsApplication::KratosPfemSolidMechanicsApplication():
    mCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mCompositeCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mCompositeCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mWallTipCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mWallTipCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mContactDomainLM2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mContactDomainPenalty2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymContactDomainLM2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymContactDomainPenalty2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )

    
  {}
  
  void KratosPfemSolidMechanicsApplication::Register()
  {
    // calling base class register to register Kratos components
    //KratosApplication::Register();
    KratosSolidMechanicsApplication::Register();
    std::cout << "Initializing KratosPfemSolidMechanicsApplication... " << std::endl;
    
    //Register Elements

    //Register Conditions

    KRATOS_REGISTER_CONDITION( "Condition2D", mCondition2D )
    KRATOS_REGISTER_CONDITION( "Condition3D", mCondition3D )

    KRATOS_REGISTER_CONDITION( "CompositeCondition2D", mCompositeCondition2D )
    KRATOS_REGISTER_CONDITION( "CompositeCondition3D", mCompositeCondition3D )
   
    KRATOS_REGISTER_CONDITION( "WallTipCondition2D", mWallTipCondition2D )
    KRATOS_REGISTER_CONDITION( "WallTipCondition3D", mWallTipCondition3D )

    KRATOS_REGISTER_CONDITION( "ContactDomainLM2DCondition", mContactDomainLM2DCondition )
    KRATOS_REGISTER_CONDITION( "ContactDomainPenalty2DCondition", mContactDomainPenalty2DCondition )

    KRATOS_REGISTER_CONDITION( "AxisymContactDomainLM2DCondition", mAxisymContactDomainLM2DCondition )
    KRATOS_REGISTER_CONDITION( "AxisymContactDomainPenalty2DCondition", mAxisymContactDomainPenalty2DCondition )

    //Register Constitutive Laws
    
    //Register Variables

    //solution
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_ACTIVE_CONTACTS )
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_STICK_CONTACTS )
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_SLIP_CONTACTS )
    
    //constitutive law	
    KRATOS_REGISTER_VARIABLE( MEAN_ERROR )

    //element

    //thermal

    //mechanical

    //geometrical
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
    KRATOS_REGISTER_VARIABLE( BOUNDARY_NORMAL )
    KRATOS_REGISTER_VARIABLE( SHRINK_FACTOR )

    //domain definition
    KRATOS_REGISTER_VARIABLE( DOMAIN_LABEL )
    KRATOS_REGISTER_VARIABLE( RIGID_WALL )
    KRATOS_REGISTER_VARIABLE( WALL_TIP_RADIUS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY )

    //contact condition
    KRATOS_REGISTER_VARIABLE( MASTER_CONDITION )
    KRATOS_REGISTER_VARIABLE( MASTER_ELEMENTS )
    KRATOS_REGISTER_VARIABLE( MASTER_NODES )

    //contact properties
    KRATOS_REGISTER_VARIABLE( FRICTION_ACTIVE )
    KRATOS_REGISTER_VARIABLE( PENALTY_PARAMETER )
    KRATOS_REGISTER_VARIABLE( TAU_STAB )
    KRATOS_REGISTER_VARIABLE( MU_STATIC )
    KRATOS_REGISTER_VARIABLE( MU_DYNAMIC )
  
  }
  
}  // namespace Kratos.


