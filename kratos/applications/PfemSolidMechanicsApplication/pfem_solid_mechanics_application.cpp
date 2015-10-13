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
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

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
  
  KRATOS_CREATE_VARIABLE(double, IMPOSED_WATER_PRESSURE )
  //constitutive law	
  KRATOS_CREATE_VARIABLE(double, MEAN_ERROR )

  //material
  KRATOS_CREATE_VARIABLE(double, PRE_CONSOLIDATION_STRESS )
  KRATOS_CREATE_VARIABLE(double, OVER_CONSOLIDATION_RATIO )
  KRATOS_CREATE_VARIABLE(double, INITIAL_SHEAR_MODULUS )
  KRATOS_CREATE_VARIABLE(double, WATER_BULK_MODULUS )
  KRATOS_CREATE_VARIABLE(double, PERMEABILITY )
  KRATOS_CREATE_VARIABLE(double, NORMAL_COMPRESSION_SLOPE )
  KRATOS_CREATE_VARIABLE(double, SWELLING_SLOPE )
  KRATOS_CREATE_VARIABLE(double, CRITICAL_STATE_LINE )
  KRATOS_CREATE_VARIABLE(double, ALPHA_SHEAR )
  KRATOS_CREATE_VARIABLE(double, INITIAL_POROSITY )
  KRATOS_CREATE_VARIABLE(double, COHESION )
  KRATOS_CREATE_VARIABLE(double, INTERNAL_DILATANCY_ANGLE )

  //element
  KRATOS_CREATE_VARIABLE(Vector, DARCY_FLOW )
  KRATOS_CREATE_VARIABLE(Matrix, TOTAL_CAUCHY_STRESS )

  // transfer variables and initial
  KRATOS_CREATE_VARIABLE(Matrix, ELASTIC_LEFT_CAUCHY_GREEN_TENSOR )
  KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_GREEN_VECTOR )

  KRATOS_CREATE_VARIABLE(Matrix, KIRCHHOFF_STRESS_TENSOR )
  KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS )

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
  KRATOS_CREATE_VARIABLE(bool,   FRICTION_ACTIVE )
  KRATOS_CREATE_VARIABLE(double, PENALTY_PARAMETER )
  KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_NORMAL )
  KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_NORMAL_REACTION )
  KRATOS_CREATE_VARIABLE(double, TAU_STAB )
  KRATOS_CREATE_VARIABLE(double, MU_STATIC )
  KRATOS_CREATE_VARIABLE(double, MU_DYNAMIC )

  // contact post process
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )
  KRATOS_CREATE_VARIABLE(double, CONTACT_ADHESION ) 
  KRATOS_CREATE_VARIABLE(double, CONTACT_FRICTION_ANGLE )
  KRATOS_CREATE_VARIABLE(double, TANGENTIAL_PENALTY_RATIO )
  KRATOS_CREATE_VARIABLE(double, CONTACT_PLASTIC_SLIP )

  // some post process variables + stress invariants
  KRATOS_CREATE_VARIABLE(double, PRECONSOLIDATION )
  KRATOS_CREATE_VARIABLE(double, STRESS_INV_P )
  KRATOS_CREATE_VARIABLE(double, STRESS_INV_J2 )
  KRATOS_CREATE_VARIABLE(double, STRESS_INV_THETA )
  KRATOS_CREATE_VARIABLE(double, VOLUMETRIC_PLASTIC )
  KRATOS_CREATE_VARIABLE(double, INCR_SHEAR_PLASTIC )

  KRATOS_CREATE_VARIABLE(double, M_MODULUS )
  KRATOS_CREATE_VARIABLE(double, SIMILAR_YOUNG_MODULUS)

  KratosPfemSolidMechanicsApplication::KratosPfemSolidMechanicsApplication():
    mTotalUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),
    mTotalUpdatedLagrangianUPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),

    // Hydro-mechanical elements
    mUpdatedLagrangianUwPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mUpdatedLagrangianUwPStabElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mUpdatedLagrangianUwPFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianUwPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianUwPStabElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),

   
    mCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mCompositeCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mCompositeCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mWallCondition2D( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mWallCondition3D( 0, Condition::GeometryType::Pointer( new Triangle3D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mContactDomainLM2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mContactDomainPenalty2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymContactDomainLM2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymContactDomainPenalty2DCondition( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )

    
  {}
  
  void KratosPfemSolidMechanicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    //KratosSolidMechanicsApplication::Register();

    std::cout << "      KRATOS' _ |  __|  __| \\   |                  " << std::endl;
    std::cout << "             '__| |_   |_   |\\ /|                  " << std::endl;
    std::cout << "            _|   _|    |__ _|  _|SOLID MECHANICS    " << std::endl;
    std::cout << "Initializing KratosPfemSolidMechanicsApplication... " << std::endl;
    
    //Register Elements

    //Register total updated lagrangian elements
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement2D3N", mTotalUpdatedLagrangianElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement2D4N", mTotalUpdatedLagrangianElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement2D6N", mTotalUpdatedLagrangianElement2D6N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement2D8N", mTotalUpdatedLagrangianElement2D8N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D4N", mTotalUpdatedLagrangianElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D6N", mTotalUpdatedLagrangianElement3D6N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D8N", mTotalUpdatedLagrangianElement3D8N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D10N", mTotalUpdatedLagrangianElement3D10N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D15N", mTotalUpdatedLagrangianElement3D15N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D20N", mTotalUpdatedLagrangianElement3D20N )
    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianElement3D27N", mTotalUpdatedLagrangianElement3D27N )

    KRATOS_REGISTER_ELEMENT( "TotalUpdatedLagrangianUPElement2D3N", mTotalUpdatedLagrangianUPElement2D3N )

    //Register updated lagrangian U wP elements
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUwPElement2D3N", mUpdatedLagrangianUwPElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUwPStabElement2D3N", mUpdatedLagrangianUwPStabElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUwPFICElement2D3N", mUpdatedLagrangianUwPFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianUwPElement2D3N", mAxisymUpdatedLagrangianUwPElement2D3N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianUwPStabElement2D3N", mAxisymUpdatedLagrangianUwPStabElement2D3N )

    //Register Conditions

    KRATOS_REGISTER_CONDITION( "Condition2D", mCondition2D )
    KRATOS_REGISTER_CONDITION( "Condition3D", mCondition3D )

    KRATOS_REGISTER_CONDITION( "CompositeCondition2D", mCompositeCondition2D )
    KRATOS_REGISTER_CONDITION( "CompositeCondition3D", mCompositeCondition3D )
   
    KRATOS_REGISTER_CONDITION( "WallCondition2D", mWallCondition2D )
    KRATOS_REGISTER_CONDITION( "WallCondition3D", mWallCondition3D )

    KRATOS_REGISTER_CONDITION( "ContactDomainLM2DCondition", mContactDomainLM2DCondition )
    KRATOS_REGISTER_CONDITION( "ContactDomainPenalty2DCondition", mContactDomainPenalty2DCondition )

    KRATOS_REGISTER_CONDITION( "AxisymContactDomainLM2DCondition", mAxisymContactDomainLM2DCondition )
    KRATOS_REGISTER_CONDITION( "AxisymContactDomainPenalty2DCondition", mAxisymContactDomainPenalty2DCondition )


    //Register Constitutive Laws
    //Serializer::Register("NonLinearHenckyCamClayPlasticPlaneStrain2DLaw", mNonLinearHenckyCamClayPlasticPlaneStrain2DLaw);
    //Serializer::Register("NonLinearHenckyCamClayPlasticAxisym2DLaw", mNonLinearHenckyCamClayPlasticAxisym2DLaw);
    //Serializer::Register("LinearHenckyCamClayPlasticPlaneStrain2DLaw", mLinearHenckyCamClayPlasticPlaneStrain2DLaw);
    //Serializer::Register("LinearHenckyCamClayPlasticAxisym2DLaw", mLinearHenckyCamClayPlasticAxisym2DLaw);
    Serializer::Register("BorjaHenckyCamClayPlasticAxisym2DLaw", mBorjaHenckyCamClayPlasticAxisym2DLaw);
    Serializer::Register("BorjaHenckyCamClayPlasticPlaneStrain2DLaw", mBorjaHenckyCamClayPlasticPlaneStrain2DLaw);
    Serializer::Register("HenckyJ2PlasticPlaneStrain2DLaw", mHenckyJ2PlasticPlaneStrain2DLaw);
    Serializer::Register("HenckyJ2PlasticAxisym2DLaw", mHenckyJ2PlasticAxisym2DLaw);
    Serializer::Register("HenckyTrescaPlasticAxisym2DLaw", mHenckyTrescaPlasticAxisym2DLaw);
    Serializer::Register("HenckyTrescaPlasticPlaneStrain2DLaw", mHenckyTrescaPlasticPlaneStrain2DLaw);
    Serializer::Register("HenckyMohrCoulombPlasticAxisym2DLaw", mHenckyMohrCoulombPlasticAxisym2DLaw);
    Serializer::Register("HenckyMohrCoulombPlasticPlaneStrain2DLaw", mHenckyMohrCoulombPlasticPlaneStrain2DLaw);

    Serializer::Register("HenckyPlasticUPJ2Axisym2DLaw", mHenckyPlasticUPJ2Axisym2DLaw);
    Serializer::Register("HenckyPlasticUPJ2PlaneStrain2DLaw", mHenckyPlasticUPJ2PlaneStrain2DLaw);
    Serializer::Register("HenckyPlasticUPTrescaAxisym2DLaw", mHenckyPlasticUPTrescaAxisym2DLaw);
    Serializer::Register("HenckyPlasticUPTrescaPlaneStrain2DLaw", mHenckyPlasticUPTrescaPlaneStrain2DLaw);

    //Register Flow Rules
    Serializer::Register("TrescaExplicitFlowRule", mTrescaExplicitFlowRule);
    Serializer::Register("J2ExplicitFlowRule", mJ2ExplicitFlowRule);
    Serializer::Register("MohrCoulombExplicitFlowRule", mMohrCoulombExplicitFlowRule);
    //Serializer::Register("CamClayExplicitFlowRule", mCamClayExplicitFlowRule);
    //Serializer::Register("LinearCamClayExplicitFlowRule", mLinearCamClayExplicitFlowRule);
    Serializer::Register("BorjaCamClayExplicitFlowRule", mBorjaCamClayExplicitFlowRule);

    //Register Yield Criterion
    Serializer::Register("J2YieldCriterion", mJ2YieldCriterion);
    Serializer::Register("TrescaYieldCriterion", mTrescaYieldCriterion);
    Serializer::Register("MohrCoulombYieldCriterion", mMohrCoulombYieldCriterion);
    Serializer::Register("CamClayYieldCriterion", mCamClayYieldCriterion);

    //Register Hardening Laws
    Serializer::Register("CamClayKinematicHardeningLaw", mCamClayKinematicHardeningLaw); 

    //Register Variables

    //solution
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_ACTIVE_CONTACTS )
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_STICK_CONTACTS )
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_SLIP_CONTACTS )


    KRATOS_REGISTER_VARIABLE( IMPOSED_WATER_PRESSURE )
    
    //constitutive law	
    KRATOS_REGISTER_VARIABLE( MEAN_ERROR )

    //material
    KRATOS_REGISTER_VARIABLE( PRE_CONSOLIDATION_STRESS )
    KRATOS_REGISTER_VARIABLE( OVER_CONSOLIDATION_RATIO )
    KRATOS_REGISTER_VARIABLE( INITIAL_SHEAR_MODULUS )
    KRATOS_REGISTER_VARIABLE( WATER_BULK_MODULUS )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY )
    KRATOS_REGISTER_VARIABLE( NORMAL_COMPRESSION_SLOPE )
    KRATOS_REGISTER_VARIABLE( SWELLING_SLOPE )
    KRATOS_REGISTER_VARIABLE( CRITICAL_STATE_LINE )
    KRATOS_REGISTER_VARIABLE( ALPHA_SHEAR )
    KRATOS_REGISTER_VARIABLE( INITIAL_POROSITY )
    KRATOS_REGISTER_VARIABLE( COHESION )
    KRATOS_REGISTER_VARIABLE( INTERNAL_DILATANCY_ANGLE ) 

    //element
    KRATOS_REGISTER_VARIABLE( DARCY_FLOW )
    KRATOS_REGISTER_VARIABLE( TOTAL_CAUCHY_STRESS )

    // transfer and initial
    KRATOS_REGISTER_VARIABLE( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS )
    KRATOS_REGISTER_VARIABLE( KIRCHHOFF_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( ELASTIC_LEFT_CAUCHY_GREEN_TENSOR )
    KRATOS_REGISTER_VARIABLE( ELASTIC_LEFT_CAUCHY_GREEN_VECTOR )

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
    KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL )
    KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL_REACTION )
    KRATOS_REGISTER_VARIABLE( TAU_STAB )
    KRATOS_REGISTER_VARIABLE( MU_STATIC )
    KRATOS_REGISTER_VARIABLE( MU_DYNAMIC )

    // Contact postprocess
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )
    KRATOS_REGISTER_VARIABLE( CONTACT_ADHESION )
    KRATOS_REGISTER_VARIABLE( CONTACT_FRICTION_ANGLE )
    KRATOS_REGISTER_VARIABLE( TANGENTIAL_PENALTY_RATIO )
    KRATOS_REGISTER_VARIABLE( CONTACT_PLASTIC_SLIP )

    // Material postprocess + invariants
    KRATOS_REGISTER_VARIABLE( PRECONSOLIDATION )
    KRATOS_REGISTER_VARIABLE( STRESS_INV_P )
    KRATOS_REGISTER_VARIABLE( STRESS_INV_J2 )
    KRATOS_REGISTER_VARIABLE( STRESS_INV_THETA )
    KRATOS_REGISTER_VARIABLE( VOLUMETRIC_PLASTIC )
    KRATOS_REGISTER_VARIABLE( INCR_SHEAR_PLASTIC )

    KRATOS_REGISTER_VARIABLE( M_MODULUS )
    KRATOS_REGISTER_VARIABLE( SIMILAR_YOUNG_MODULUS )

  }
  
}  // namespace Kratos.


