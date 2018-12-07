//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:               January 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/mat_variables.h"
#include "custom_friction/friction_law.hpp"

#include "custom_conditions/hydraulic_contact/hydraulic_axisym_rigid_contact_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

   //************************************************************************************
   //************************************************************************************

   HydraulicAxisymRigidContactPenalty2DCondition::HydraulicAxisymRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
      : HydraulicRigidContactPenalty3DCondition(NewId, pGeometry)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   HydraulicAxisymRigidContactPenalty2DCondition::HydraulicAxisymRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : HydraulicRigidContactPenalty3DCondition(NewId, pGeometry, pProperties)
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //************************************************************************************
   //************************************************************************************

   HydraulicAxisymRigidContactPenalty2DCondition::HydraulicAxisymRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : HydraulicRigidContactPenalty3DCondition(NewId, pGeometry, pProperties, pRigidWall)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   HydraulicAxisymRigidContactPenalty2DCondition::HydraulicAxisymRigidContactPenalty2DCondition( HydraulicAxisymRigidContactPenalty2DCondition const& rOther )
      : HydraulicRigidContactPenalty3DCondition(rOther)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer HydraulicAxisymRigidContactPenalty2DCondition::Create(IndexType NewId, const NodesArrayType& ThisNodes, PropertiesType::Pointer pProperties) const
   {
     return Kratos::make_shared<HydraulicAxisymRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
   }

   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer HydraulicAxisymRigidContactPenalty2DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
   {
     return Kratos::make_shared<HydraulicAxisymRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
   }


   //************************************************************************************
   //************************************************************************************

   HydraulicAxisymRigidContactPenalty2DCondition::~HydraulicAxisymRigidContactPenalty2DCondition()
   {

   }



   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************

   void HydraulicAxisymRigidContactPenalty2DCondition::CalculateKinematics(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo, const double& rPointNumber)
   {
      KRATOS_TRY

      CalculateRadius( rVariables.CurrentRadius, rVariables.ReferenceRadius );

      HydraulicRigidContactPenalty3DCondition::CalculateKinematics( rVariables, rCurrentProcessInfo, rPointNumber);

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void  HydraulicAxisymRigidContactPenalty2DCondition::CalculateRadius(double & rCurrentRadius, double & rReferenceRadius)
   {

      KRATOS_TRY

      rCurrentRadius=0;
      rReferenceRadius=0;

      //Displacement from the reference to the current configuration
      array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
      array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);
      array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
      array_1d<double, 3 > & CurrentPosition      = GetGeometry()[0].Coordinates();
      array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;


      rCurrentRadius   = CurrentPosition[0];
      rReferenceRadius = ReferencePosition[0];

      if( rCurrentRadius == 0 ){

         WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

         double counter = 0;

         for(unsigned int i = 0; i < rN.size(); i++)
         {
            array_1d<double, 3 > & NodePosition = rN[i].Coordinates();
            if( NodePosition[0] != 0 ){
               rCurrentRadius += NodePosition[0] * 0.225;
               counter ++;
            }

         }

         rCurrentRadius /= counter;

      }


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void  HydraulicAxisymRigidContactPenalty2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
   {

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

      if ( this->GetProperties().Has(THICKNESS) ) {
         if( GetProperties()[THICKNESS] > 0 )
            IntegrationWeight /=  GetProperties()[THICKNESS];
      }

      //contributions to stiffness matrix calculated on the reference config

      HydraulicRigidContactPenalty3DCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

      //KRATOS_WATCH( rLeftHandSideMatrix )
   }


   //************************************************************************************
   //************************************************************************************

   void   HydraulicAxisymRigidContactPenalty2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
   {
      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

      if ( this->GetProperties().Has(THICKNESS) ) {
         if( GetProperties()[THICKNESS] > 0 )
            IntegrationWeight /=  GetProperties()[THICKNESS];
      }

      //contribution to external forces

      HydraulicRigidContactPenalty3DCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

      //KRATOS_WATCH( rRightHandSideVector )
   }






} // Namespace Kratos
