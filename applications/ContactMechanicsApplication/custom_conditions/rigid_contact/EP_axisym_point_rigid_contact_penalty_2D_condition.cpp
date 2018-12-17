//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/rigid_contact/EP_axisym_point_rigid_contact_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{
   //************************************************************************************
   //************************************************************************************
   EPAxisymPointRigidContactPenalty2DCondition::EPAxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer
         pGeometry)
      : EPPointRigidContactPenalty2DCondition(NewId, pGeometry)
   {
      //DO NOT ADD DOFS HERE!!!

   }

   //************************************************************************************
   //************************************************************************************
   EPAxisymPointRigidContactPenalty2DCondition::EPAxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : EPPointRigidContactPenalty2DCondition(NewId, pGeometry, pProperties)
   {
   }


   //************************************************************************************
   //************************************************************************************
   EPAxisymPointRigidContactPenalty2DCondition::EPAxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : EPPointRigidContactPenalty2DCondition(NewId, pGeometry, pProperties, pRigidWall)
   {

   }

   //************************************************************************************
   //************************************************************************************
   EPAxisymPointRigidContactPenalty2DCondition::EPAxisymPointRigidContactPenalty2DCondition( EPAxisymPointRigidContactPenalty2DCondition const& rOther )
      : EPPointRigidContactPenalty2DCondition(rOther)
   {
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer EPAxisymPointRigidContactPenalty2DCondition::Create(IndexType NewId, NodesArrayType
         const& ThisNodes,  PropertiesType::Pointer pProperties) const
   {
     return Kratos::make_shared<EPAxisymPointRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
   }

   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer EPAxisymPointRigidContactPenalty2DCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
   {
      EPAxisymPointRigidContactPenalty2DCondition NewCondition( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
      NewCondition.mCurrentInfo = this->mCurrentInfo;
      NewCondition.mSavedInfo   = this->mSavedInfo;

      // in the constructor of NewCondition I create a new friction law and here I clone the this->
      NewCondition.mpFrictionLaw = this->mpFrictionLaw->Clone();

      return Kratos::make_shared<EPAxisymPointRigidContactPenalty2DCondition>(NewCondition);
   }

   //************************************************************************************
   //************************************************************************************


   EPAxisymPointRigidContactPenalty2DCondition::~EPAxisymPointRigidContactPenalty2DCondition()
   {
   }



   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************

   void EPAxisymPointRigidContactPenalty2DCondition::CalculateKinematics(ConditionVariables& rVariables,
         const ProcessInfo& rCurrentProcessInfo,
         const double& rPointNumber)
   {
      KRATOS_TRY

      EPPointRigidContactPenalty2DCondition::CalculateKinematics(rVariables, rCurrentProcessInfo, rPointNumber);

      CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius);

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void  EPAxisymPointRigidContactPenalty2DCondition::CalculateRadius(double & rCurrentRadius, double & rReferenceRadius)
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

   void  EPAxisymPointRigidContactPenalty2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
   {

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

      if ( this->GetProperties().Has(THICKNESS) ) {
         if( GetProperties()[THICKNESS] > 0 )
            IntegrationWeight /=  GetProperties()[THICKNESS];
      }

      //contributions to stiffness matrix calculated on the reference config

      EPPointRigidContactPenalty2DCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

      //KRATOS_WATCH( rLeftHandSideMatrix )
   }


   //************************************************************************************
   //************************************************************************************

   void  EPAxisymPointRigidContactPenalty2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
   {
      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

      if ( this->GetProperties().Has(THICKNESS) ) {
         if( GetProperties()[THICKNESS] > 0 )
            IntegrationWeight /=  GetProperties()[THICKNESS];
      }

      //contribution to external forces

      EPPointRigidContactPenalty2DCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

      //KRATOS_WATCH( rRightHandSideVector )
   }


} // Namespace Kratos
