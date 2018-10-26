//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

#include "contact_mechanics_application.h"

namespace Kratos {

  //Application variables creation: (see contact_mechanics_application_variables.cpp)

  //Application Constructor:
  KratosContactMechanicsApplication::KratosContactMechanicsApplication():
    KratosApplication("ContactMechanicsApplication"),
    mContactDomainLMCondition3D4N( 0, Kratos::make_shared< Tetrahedra3D4<Node<3> > >( Condition::GeometryType::PointsArrayType(4))),
    mContactDomainLMCondition2D3N( 0, Kratos::make_shared< Triangle2D3<Node<3> > >( Condition::GeometryType::PointsArrayType(3))),
    mContactDomainPenaltyCondition2D3N( 0, Kratos::make_shared< Triangle2D3<Node<3> > >( Condition::GeometryType::PointsArrayType(3))),
    mAxisymContactDomainLMCondition2D3N( 0, Kratos::make_shared< Triangle2D3<Node<3> > >( Condition::GeometryType::PointsArrayType(3))),
    mAxisymContactDomainPenaltyCondition2D3N( 0, Kratos::make_shared< Triangle2D3<Node<3> > >( Condition::GeometryType::PointsArrayType(3))),
    mThermalContactDomainPenaltyCondition2D3N( 0, Kratos::make_shared< Triangle2D3<Node<3> > >( Condition::GeometryType::PointsArrayType(3))),
    mAxisymThermalContactDomainPenaltyCondition2D3N( 0, Kratos::make_shared< Triangle2D3<Node<3> > >( Condition::GeometryType::PointsArrayType(3)))
    // ,
    // mRigidBodyPointLinkCondition2D1N( 0, Kratos::make_shared< Point3D<Node<3> >( Condition::GeometryType::PointsArrayType(1))),
    // mRigidBodyPointLinkCondition3D1N( 0, Kratos::make_shared< Point3D<Node<3> >( Condition::GeometryType::PointsArrayType(1)))

  {}

  void KratosContactMechanicsApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::cout << "             ___         _           _          " << std::endl;
    std::cout << "     KRATOS / __|___ _ _| |_ __ _ __| |_          " << std::endl;
    std::cout << "           | (__/ _ \\ ' \\  _/ _` / _|  _|         " << std::endl;
    std::cout << "            \\___\\___/_||_\\__\\__,_\\__|\\__|MECHANICS" << std::endl;
    std::cout << "Initializing KratosContactMechanicsApplication... " << std::endl;


      //Register Rigid Bodies
      // Serializer::Register( "RigidBodyElement", mRigidBodyElement);
      // Serializer::Register( "TranslatoryRigidBodyElement", mTranslatoryRigidBodyElement);
      KRATOS_REGISTER_ELEMENT( "RigidBodyElement", mRigidBodyElement )
      KRATOS_REGISTER_ELEMENT( "TranslatoryRigidBodyElement", mTranslatoryRigidBodyElement )
      KRATOS_REGISTER_ELEMENT( "RigidBodySegregatedVElement", mRigidBodySegregatedVElement )

      //Register Conditions
      KRATOS_REGISTER_CONDITION( "ContactDomainLMCondition3D4N", mContactDomainLMCondition3D4N )

      KRATOS_REGISTER_CONDITION( "ContactDomainLMCondition2D3N", mContactDomainLMCondition2D3N )
      KRATOS_REGISTER_CONDITION( "ContactDomainPenaltyCondition2D3N", mContactDomainPenaltyCondition2D3N )

      KRATOS_REGISTER_CONDITION( "AxisymContactDomainLMCondition2D3N", mAxisymContactDomainLMCondition2D3N )
      KRATOS_REGISTER_CONDITION( "AxisymContactDomainPenaltyCondition2D3N", mAxisymContactDomainPenaltyCondition2D3N )

      KRATOS_REGISTER_CONDITION( "ThermalContactDomainPenaltyCondition2D3N", mThermalContactDomainPenaltyCondition2D3N )
      KRATOS_REGISTER_CONDITION( "AxisymThermalContactDomainPenaltyCondition2D3N", mAxisymThermalContactDomainPenaltyCondition2D3N )

      // KRATOS_REGISTER_CONDITION( "RigidBodyPointLinkCondition2D1N", mRigidBodyPointLinkCondition2D1N )
      // KRATOS_REGISTER_CONDITION( "RigidBodyPointLinkCondition3D1N", mRigidBodyPointLinkCondition3D1N )

      // Serializer::Register( "PointRigidContactPenalty2DCondition", mPointRigidContactPenalty2DCondition );
      // Serializer::Register( "PointRigidContactPenalty3DCondition", mPointRigidContactPenalty3DCondition );
      // Serializer::Register( "AxisymPointRigidContactPenalty2DCondition", mAxisymPointRigidContactPenalty2DCondition );

      // Serializer::Register( "EPPointRigidContactPenalty2DCondition", mEPPointRigidContactPenalty2DCondition );
      // Serializer::Register( "EPPointRigidContactPenalty3DCondition", mEPPointRigidContactPenalty3DCondition );
      // Serializer::Register( "EPAxisymPointRigidContactPenalty2DCondition", mEPAxisymPointRigidContactPenalty2DCondition );

      // Serializer::Register( "HydraulicRigidContactPenalty3DCondition", mHydraulicRigidContactPenalty3DCondition );
      // Serializer::Register( "HydraulicAxisymRigidContactPenalty2DCondition", mHydraulicAxisymRigidContactPenalty2DCondition );
      KRATOS_REGISTER_CONDITION( "PointRigidContactPenalty2DCondition", mPointRigidContactPenalty2DCondition );
      KRATOS_REGISTER_CONDITION( "PointRigidContactPenalty3DCondition", mPointRigidContactPenalty3DCondition );
      KRATOS_REGISTER_CONDITION( "AxisymPointRigidContactPenalty2DCondition", mAxisymPointRigidContactPenalty2DCondition );

      KRATOS_REGISTER_CONDITION( "EPPointRigidContactPenalty2DCondition", mEPPointRigidContactPenalty2DCondition );
      KRATOS_REGISTER_CONDITION( "EPPointRigidContactPenalty3DCondition", mEPPointRigidContactPenalty3DCondition );
      KRATOS_REGISTER_CONDITION( "EPAxisymPointRigidContactPenalty2DCondition", mEPAxisymPointRigidContactPenalty2DCondition );

      KRATOS_REGISTER_CONDITION( "HydraulicRigidContactPenalty3DCondition", mHydraulicRigidContactPenalty3DCondition );
      KRATOS_REGISTER_CONDITION( "HydraulicAxisymRigidContactPenalty2DCondition", mHydraulicAxisymRigidContactPenalty2DCondition );


      //Register friction laws
      Serializer::Register( "FrictionLaw", mFrictionLaw );
      Serializer::Register( "CoulombAdhesionFrictionLaw", mCoulombAdhesionFrictionLaw );
      Serializer::Register( "HardeningCoulombFrictionLaw", mHardeningCoulombFrictionLaw );

      //Register Variables
      KRATOS_REGISTER_VARIABLE( FRICTION_LAW_NAME )
      KRATOS_REGISTER_VARIABLE( FRICTION_LAW )
      KRATOS_REGISTER_VARIABLE( HYDRAULIC )

      //contact properties
      KRATOS_REGISTER_VARIABLE( FRICTION_ACTIVE )
      KRATOS_REGISTER_VARIABLE( PENALTY_PARAMETER )
      KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL )
      KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL_REACTION )
      KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_TANGENTIAL )
      KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_TANGENTIAL_REACTION )
      KRATOS_REGISTER_VARIABLE( TAU_STAB )
      KRATOS_REGISTER_VARIABLE( MU_STATIC )
      KRATOS_REGISTER_VARIABLE( MU_DYNAMIC )

      //contact postprocess
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WATER_CONTACT_FORCE )
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )
      KRATOS_REGISTER_VARIABLE( CONTACT_ADHESION )
      KRATOS_REGISTER_VARIABLE( CONTACT_FRICTION_ANGLE )
      KRATOS_REGISTER_VARIABLE( TANGENTIAL_PENALTY_RATIO )
      KRATOS_REGISTER_VARIABLE( CONTACT_PLASTIC_SLIP )

      //thermal properties
      KRATOS_REGISTER_VARIABLE( HEAT_CONDUCTIVITY )

      //solution
      KRATOS_REGISTER_VARIABLE(SEGREGATED_STEP)    
      KRATOS_REGISTER_VARIABLE(CONTACT_STEP_TIME)

      }

}  // namespace Kratos.
