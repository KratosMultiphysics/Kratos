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
#include "custom_conditions/hydraulic_contact/hydraulic_rigid_contact_penalty_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

   //************************************************************************************
   //************************************************************************************

   HydraulicRigidContactPenalty3DCondition::HydraulicRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
      : PointRigidContactCondition(NewId, pGeometry)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   HydraulicRigidContactPenalty3DCondition::HydraulicRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : PointRigidContactCondition(NewId, pGeometry, pProperties)
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //************************************************************************************
   //************************************************************************************

   HydraulicRigidContactPenalty3DCondition::HydraulicRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : PointRigidContactCondition(NewId, pGeometry, pProperties, pRigidWall)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   HydraulicRigidContactPenalty3DCondition::HydraulicRigidContactPenalty3DCondition( HydraulicRigidContactPenalty3DCondition const& rOther )
      : PointRigidContactCondition(rOther)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer HydraulicRigidContactPenalty3DCondition::Create(IndexType NewId, const NodesArrayType& ThisNodes, PropertiesType::Pointer pProperties) const
   {
     return Kratos::make_shared<HydraulicRigidContactPenalty3DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
   }

   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer HydraulicRigidContactPenalty3DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
   {
     return Kratos::make_shared<HydraulicRigidContactPenalty3DCondition>(NewId,GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
   }


   //************************************************************************************
   //************************************************************************************

   HydraulicRigidContactPenalty3DCondition::~HydraulicRigidContactPenalty3DCondition()
   {

   }

   void HydraulicRigidContactPenalty3DCondition::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      // virtual because I don't want to delete "mechanical" contact forces
      if ( GetGeometry()[0].SolutionStepsDataHas( WATER_CONTACT_FORCE )) {
         array_1d<double, 3 > & rWaterContactForce = GetGeometry()[0].FastGetSolutionStepValue( WATER_CONTACT_FORCE );

         for(unsigned int j = 0; j < 3; j++)
         {
            rWaterContactForce[j] = 0;
         }
      }

      KRATOS_CATCH("")
   }

   //************* GETTING METHODS

   //***********************************************************************************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::GetDofList(DofsVectorType& rConditionDofList,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      rConditionDofList.resize(0);
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(WATER_DISPLACEMENT_X));
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(WATER_DISPLACEMENT_Y));
         if( dimension == 3 )
            rConditionDofList.push_back(GetGeometry()[i].pGetDof(WATER_DISPLACEMENT_Z));
      }


      KRATOS_CATCH( "" )
   }

   //***********************************************************************************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::EquationIdVector(EquationIdVectorType& rResult,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int condition_size        = number_of_nodes * dimension;

      if (rResult.size() != condition_size)
         rResult.resize( condition_size, false );

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         int index = i * dimension;
         rResult[index]     = GetGeometry()[i].GetDof(WATER_DISPLACEMENT_X).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof(WATER_DISPLACEMENT_Y).EquationId();
         if( dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof(WATER_DISPLACEMENT_Z).EquationId();
      }

      KRATOS_CATCH( "" )
   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::GetValuesVector(Vector& rValues, int Step)
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size  = number_of_nodes * dimension;

      if ( rValues.size() != condition_size )
         rValues.resize( condition_size, false );

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         unsigned int index = i * dimension;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Y, Step );

         if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Z, Step );
      }
   }

   //***********************************************************************************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size    = number_of_nodes * dimension;

      if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Y, Step );

         if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Z, Step );
      }
   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size    = number_of_nodes * dimension;

      if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Y, Step );

         if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Z, Step );
      }

   }

   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************

   void HydraulicRigidContactPenalty3DCondition::CalculateKinematics(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo, const double& rPointNumber)
   {
      KRATOS_TRY

      SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, rVariables.RelativeDisplacement);


      if( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ){

         rVariables.Options.Set(ACTIVE,true);

         rVariables.Gap.Normal = fabs(rVariables.Gap.Normal);

         //get contact properties and parameters
         this->CalculateContactFactors( rVariables );

      }
      else{

         rVariables.Options.Set(ACTIVE,false);
      }

      rVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void HydraulicRigidContactPenalty3DCondition::CalculateContactFactors(ConditionVariables &rVariables)
   {

      KRATOS_TRY

      //Compute the neighbour distance, then a stress-"like" may be computed.
      WeakPointerVector<Node<3> >& rN  = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);
      array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
      array_1d<double,3> Neighb_Point;

      double distance = 0;
      double counter = 0;

      for(unsigned int i = 0; i < rN.size(); i++)
      {
         if(rN[i].Is(BOUNDARY)){

            Neighb_Point[0] = rN[i].X();
            Neighb_Point[1] = rN[i].Y();
            Neighb_Point[2] = rN[i].Z();

            distance += norm_2(Contact_Point-Neighb_Point);

            counter ++;
         }
      }

      if( counter != 0 )
         distance /= counter;

      if( distance == 0 )
         distance = 1;


      rVariables.ContributoryFactor = distance;

      //get contact properties and parameters
      double PenaltyParameter = 1;
      if( GetProperties().Has(PENALTY_PARAMETER) )
         PenaltyParameter = GetProperties()[PENALTY_PARAMETER];

      WeakPointerVector<Element >& rE = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
      double ElasticModulus = 0;
      if( GetProperties().Has(YOUNG_MODULUS) )
         ElasticModulus = GetProperties()[YOUNG_MODULUS];
      else
         ElasticModulus = rE.front().GetProperties()[YOUNG_MODULUS];

      // the Modified Cam Clay model does not have a constant Young modulus, so something similar to that is computed
      if (ElasticModulus <= 1.0e-5) {
         std::vector<double> mModulus;
         ProcessInfo SomeProcessInfo;
         for ( unsigned int i = 0; i < rE.size(); i++)
         {
            rE[i].CalculateOnIntegrationPoints(EQUIVALENT_YOUNG_MODULUS, mModulus, SomeProcessInfo);
            ElasticModulus += mModulus[0];
         }
         ElasticModulus /= double(rE.size());
      }

      double factor = 4;
      if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
         int order = (int)((-1) * std::log10(distance) + 1) ;
         distance *= factor * pow(10,order);
      }

      rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;

      //std::cout<<" ContactPoint["<<this->Id()<<"]: penalty_n"<<rVariables.Penalty.Normal<<", ElasticModulus: "<<ElasticModulus<<", distance: "<<distance<<std::endl;

      KRATOS_CATCH( "" )

   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      if( rVariables.Options.Is(ACTIVE)){

         Matrix Kuug(3,3);
         noalias(Kuug) = ZeroMatrix(3,3);

         noalias(Kuug) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);

         for(unsigned int i=0; i<dimension; i++)
         {
            for(unsigned int j=0; j<dimension; j++)
            {
               rLeftHandSideMatrix(i,j) = Kuug(i,j);
            }
         }


      }
      else{

         rLeftHandSideMatrix= ZeroMatrix(dimension,dimension);
      }

      // if( rVariables.Options.Is(ACTIVE))
      //   std::cout<<" Contact Tangent Matrix ["<<this->Id()<<"]: "<<rLeftHandSideMatrix<<std::endl;

      KRATOS_CATCH( "" )
   }


   void HydraulicRigidContactPenalty3DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      if( rVariables.Options.Is(ACTIVE)){

         this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

      }
      else{

         rRightHandSideVector = ZeroVector(dimension);

      }

      // if( rVariables.Options.Is(ACTIVE)){
      //   std::cout<<" Contact Forces Vector ["<<this->Id()<<"]: "<<rRightHandSideVector<<std::endl;
      //   std::cout<<" Tangent Force "<<GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE)<<std::endl;

      // }

      KRATOS_CATCH( "" )
   }

   //**************************** Calculate Normal Contact Force ***********************
   //***********************************************************************************

   void HydraulicRigidContactPenalty3DCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {

      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

      NormalForceModulus *= rIntegrationWeight;


      for(unsigned int j = 0; j < dimension; j++)
      {
         rRightHandSideVector[j] = -NormalForceModulus * rVariables.Surface.Normal[j];
      }

      if ( GetGeometry()[0].SolutionStepsDataHas( WATER_CONTACT_FORCE )) {
         array_1d<double, 3 > & rWaterContactForce = GetGeometry()[0].FastGetSolutionStepValue( WATER_CONTACT_FORCE );
         for(unsigned int j = 0; j < dimension; j++)
         {
            rWaterContactForce[j] = NormalForceModulus * rVariables.Surface.Normal[j];
         }
      }

      KRATOS_CATCH("")

   }


   //**************************** Calculate Normal Force Modulus ***********************
   //***********************************************************************************

   double& HydraulicRigidContactPenalty3DCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, ConditionVariables& rVariables )
   {
      KRATOS_TRY


      const array_1d< double, 3> & CurrentWaterDisplacement = GetGeometry()[0].FastGetSolutionStepValue( WATER_DISPLACEMENT);
      const array_1d< double, 3> & PreviousWaterDisplacement = GetGeometry()[0].FastGetSolutionStepValue( WATER_DISPLACEMENT, 1);

      array_1d<double, 3> DeltaWaterPressure;
      DeltaWaterPressure = CurrentWaterDisplacement - PreviousWaterDisplacement;

      rNormalForceModulus = 0;
      for (unsigned int i = 0; i < 3; i++) {
         rNormalForceModulus += DeltaWaterPressure[i] * rVariables.Surface.Normal(i);
      }

      rNormalForceModulus *= rVariables.Penalty.Normal;

      return rNormalForceModulus;

      KRATOS_CATCH( "" )

   }




} // Namespace Kratos
