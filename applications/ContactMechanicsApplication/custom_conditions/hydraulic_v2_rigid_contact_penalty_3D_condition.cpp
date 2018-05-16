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
#include "custom_conditions/hydraulic_v2_rigid_contact_penalty_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"

#include "custom_friction/friction_law.hpp"

#include "includes/mat_variables.h"


namespace Kratos
{

   //************************************************************************************
   //************************************************************************************

   HydraulicV2RigidContactPenalty3DCondition::HydraulicV2RigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
      : HydraulicRigidContactPenalty3DCondition(NewId, pGeometry)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   HydraulicV2RigidContactPenalty3DCondition::HydraulicV2RigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : HydraulicRigidContactPenalty3DCondition(NewId, pGeometry, pProperties)
   {
      //DO NOT ADD DOFS HERE!!!    
   }


   //************************************************************************************
   //************************************************************************************

   HydraulicV2RigidContactPenalty3DCondition::HydraulicV2RigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : HydraulicRigidContactPenalty3DCondition(NewId, pGeometry, pProperties, pRigidWall)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   HydraulicV2RigidContactPenalty3DCondition::HydraulicV2RigidContactPenalty3DCondition( HydraulicV2RigidContactPenalty3DCondition const& rOther )
      : HydraulicRigidContactPenalty3DCondition(rOther)
   {
      //DO NOT ADD DOFS HERE!!! 
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer HydraulicV2RigidContactPenalty3DCondition::Create(IndexType NewId, const NodesArrayType& ThisNodes, PropertiesType::Pointer pProperties) const
   {
      return Condition::Pointer(new HydraulicV2RigidContactPenalty3DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
   }

   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer HydraulicV2RigidContactPenalty3DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
   {
      return Condition::Pointer(new HydraulicV2RigidContactPenalty3DCondition(NewId,GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall));
   }


   //************************************************************************************
   //************************************************************************************

   HydraulicV2RigidContactPenalty3DCondition::~HydraulicV2RigidContactPenalty3DCondition()
   {

   }

   //************* GETTING METHODS

   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::GetDofList(DofsVectorType& rConditionDofList,
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
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(WATER_PRESSURE));
      }


      KRATOS_CATCH( "" )
   }

   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::EquationIdVector(EquationIdVectorType& rResult,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int condition_size        = number_of_nodes * (dimension + 1);

      if (rResult.size() != condition_size)
         rResult.resize( condition_size, false );

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         int index = i * dimension;
         rResult[index]     = GetGeometry()[i].GetDof(WATER_DISPLACEMENT_X).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof(WATER_DISPLACEMENT_Y).EquationId();
         if( dimension == 3) {
            rResult[index + 2] = GetGeometry()[i].GetDof(WATER_DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
         }
         rResult[index + 2] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
      }

      KRATOS_CATCH( "" )
   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::GetValuesVector(Vector& rValues, int Step)
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size  = number_of_nodes * ( dimension + 1);

      if ( rValues.size() != condition_size ) 
         rValues.resize( condition_size, false );

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         unsigned int index = i * dimension;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Y, Step );

         if ( dimension == 3 ) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
      }
   }

   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size    = number_of_nodes * (dimension + 1);

      if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Y, Step );

         if ( dimension == 3 ) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Z, Step );
            rValues[index + 3] = 0.0;
         }
         rValues[index + 2] = 0.0;
      }
   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size    = number_of_nodes * (dimension + 1);

      if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Y, Step );

         if ( dimension == 3 ) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Z, Step );
            rValues[index + 3] = 0.0;
         }
         rValues[index + 2] = 0.0;
      }

   }

   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::CalculateKinematics(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo, const double& rPointNumber)
   {
      KRATOS_TRY

      SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, rVariables.RelativeDisplacement);


      mDirichletCondition = false;


      if( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ){

         rVariables.Options.Set(ACTIVE,true);

         rVariables.Gap.Normal = fabs(rVariables.Gap.Normal);

         //get contact properties and parameters
         this->CalculateContactFactors( rVariables );

      }
      else{

         rVariables.Options.Set(ACTIVE,false);
         const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();


         mDirichletCondition = true;
        
         if ( GetGeometry()[0].pGetDof(WATER_DISPLACEMENT_X)->IsFixed() )
            mDirichletCondition = false;
         if ( GetGeometry()[0].pGetDof(WATER_DISPLACEMENT_Y)->IsFixed() ) {
            mDirichletCondition = false;
         }
         if ( GetGeometry()[0].HasDofFor(WATER_DISPLACEMENT_Z) ) {
            if ( GetGeometry()[0].pGetDof(WATER_DISPLACEMENT_Z)->IsFixed() && dimension > 2) {
               mDirichletCondition = false;
            }
         }
         if ( GetGeometry()[0].pGetDof(WATER_PRESSURE)->IsFixed() )
            mDirichletCondition = false;

         if (mDirichletCondition) {

            //get contact properties and parameters
            this->CalculateContactFactors( rVariables );

            rVariables.Penalty.Normal = 1.0e4;
         }

      }

      rVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * (dimension + 1);

      if ( rCalculationFlags.Is(HydraulicRigidContactPenalty3DCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
         if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

         noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }


      //resizing as needed the RHS
      if ( rCalculationFlags.Is(HydraulicRigidContactPenalty3DCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
         if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

         rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

      }
   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
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
      else if ( mDirichletCondition) {

         // DO SOMETHING
         rLeftHandSideMatrix(dimension, dimension) =  rVariables.Penalty.Normal *  rIntegrationWeight;
      }
      else{

         rLeftHandSideMatrix= ZeroMatrix(dimension+1,dimension+1);
      }

      // if( rVariables.Options.Is(ACTIVE))
      //   std::cout<<" Contact Tangent Matrix ["<<this->Id()<<"]: "<<rLeftHandSideMatrix<<std::endl;

      KRATOS_CATCH( "" )
   }


   //***********************************************************************************
   //***********************************************************************************

   void HydraulicV2RigidContactPenalty3DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      if( rVariables.Options.Is(ACTIVE)){

         this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

      }
      else if ( mDirichletCondition) {
         // DO SOMETHING
         const double & rWaterPressure = GetGeometry()[0].FastGetSolutionStepValue( WATER_PRESSURE);
         rRightHandSideVector[dimension] = -rVariables.Penalty.Normal * rIntegrationWeight * rWaterPressure;
         //std::cout << " this->Id() " << this->Id() << " , " << rVariables.Penalty.Normal << " , " << rIntegrationWeight << " , " << rWaterPressure << std::endl;
         //std::cout << " " <<  rRightHandSideVector << std::endl;
      }
      else{

         rRightHandSideVector = ZeroVector(dimension+1);

      }

      // if( rVariables.Options.Is(ACTIVE)){
      //   std::cout<<" Contact Forces Vector ["<<this->Id()<<"]: "<<rRightHandSideVector<<std::endl;
      //   std::cout<<" Tangent Force "<<GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE)<<std::endl;

      // }

      KRATOS_CATCH( "" )
   }



} // Namespace Kratos



