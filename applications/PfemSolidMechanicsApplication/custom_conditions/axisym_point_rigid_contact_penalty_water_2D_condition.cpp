//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//


#include "axisym_point_rigid_contact_penalty_water_2D_condition.hpp"

#include "pfem_solid_mechanics_application_variables.h"

// This is the same than  AxisymPointRigidContactPenalty2DCondition; however, if the yield surface of the contact is written with in terms of the 
// effective contact stress, it induces a new term in the global matrix. (THIS IS WHY IT IS PROGRAMED; just to check if it increases the convergence 
// of the problem

namespace Kratos
{
   // USUAL CONSTRUCTORS

   //************************************************************************************
   //************************************************************************************
   AxisymPointRigidContactPenaltyWater2DCondition::AxisymPointRigidContactPenaltyWater2DCondition(IndexType NewId, GeometryType::Pointer
         pGeometry)
      : AxisymPointRigidContactPenalty2DCondition(NewId, pGeometry)
   {
   }

   //************************************************************************************
   //************************************************************************************
   AxisymPointRigidContactPenaltyWater2DCondition::AxisymPointRigidContactPenaltyWater2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : AxisymPointRigidContactPenalty2DCondition(NewId, pGeometry, pProperties)
   {
   }


   //************************************************************************************
   //************************************************************************************
   AxisymPointRigidContactPenaltyWater2DCondition::AxisymPointRigidContactPenaltyWater2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : AxisymPointRigidContactPenalty2DCondition(NewId, pGeometry, pProperties, pRigidWall)
   {

   }

   //************************************************************************************
   //************************************************************************************
   AxisymPointRigidContactPenaltyWater2DCondition::AxisymPointRigidContactPenaltyWater2DCondition( AxisymPointRigidContactPenaltyWater2DCondition const& rOther )
      : AxisymPointRigidContactPenalty2DCondition(rOther)
   {
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer AxisymPointRigidContactPenaltyWater2DCondition::Create(IndexType NewId, NodesArrayType
         const& ThisNodes,  PropertiesType::Pointer pProperties) const
   {
      return Condition::Pointer(new AxisymPointRigidContactPenaltyWater2DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
   }


   //************************************************************************************
   //************************************************************************************


   AxisymPointRigidContactPenaltyWater2DCondition::~AxisymPointRigidContactPenaltyWater2DCondition()
   {
   }

   // AQUI VIENEN LAS COSAS

   // ADD A NEW DOF DUE TO THE WATER
   void AxisymPointRigidContactPenaltyWater2DCondition::GetDofList(DofsVectorType& rConditionDofList,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      rConditionDofList.resize(0);
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(WATER_PRESSURE));
      }


      KRATOS_CATCH( "" )
   }

   //***********************************************************************************
   //***********************************************************************************

   void AxisymPointRigidContactPenaltyWater2DCondition::EquationIdVector(EquationIdVectorType& rResult,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int condition_size        = number_of_nodes * (dimension + 1); // due to the water pressure

      if (rResult.size() != condition_size)
         rResult.resize( condition_size, false );

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         int index = i * (dimension + 1);
         rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
         rResult[index + 2] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
      }

      KRATOS_CATCH( "" )
   }


   //***********************************************************************************
   //***********************************************************************************

   void AxisymPointRigidContactPenaltyWater2DCondition::GetValuesVector(Vector& rValues, int Step)
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int condition_size        = number_of_nodes * (dimension + 1); // due to the water pressure

      if ( rValues.size() != condition_size ) 
         rValues.resize( condition_size, false );

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         unsigned int index = i * (dimension+1);
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
         rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
      }
   }

   //***********************************************************************************
   //***********************************************************************************

   void AxisymPointRigidContactPenaltyWater2DCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int condition_size        = number_of_nodes * (dimension + 1); // due to the water pressure

      if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension + 1);
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
         rValues[index + 2] = 0.0;

      }
   }


   //***********************************************************************************
   //***********************************************************************************

   void AxisymPointRigidContactPenaltyWater2DCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int condition_size        = number_of_nodes * (dimension + 1); // due to the water pressure

      if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension + 1);
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
         rValues[index + 2] = 0.0;
      }

   }


   // MODIFY THE INIZIALIZATION due to the differnt size of the system due to the water pressure
   void AxisymPointRigidContactPenaltyWater2DCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

      if ( rCalculationFlags.Is(PointRigidContactCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
         if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

         noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }


      //resizing as needed the RHS
      if ( rCalculationFlags.Is(PointRigidContactCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
         if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

         rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

      }
   }


   // COMPUTE Kuug. ( the same, but since is a different size, is copied and put it ok)
   void AxisymPointRigidContactPenaltyWater2DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {

      KRATOS_TRY
      // it is better to re-copy the normal part, and then re-do the tangential part with the water-pressure part....
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;
      rLeftHandSideMatrix= ZeroMatrix(MatSize,MatSize);   

      if ( rVariables.Options.Is(ACTIVE) )
      {
         MatrixType  LocalKuugNormal = ZeroMatrix(MatSize-1, MatSize-1);
         LocalKuugNormal = rVariables.Penalty.Normal* rIntegrationWeight * outer_prod( rVariables.Surface.Normal, rVariables.Surface.Normal);
         for (unsigned int i = 0; i < MatSize-1; i++){
            for (unsigned int j = 0; j < MatSize-1; j++){
               rLeftHandSideMatrix(i,j) = LocalKuugNormal(i,j);
            }
         }
         this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight);  // MINE THAT HAS THE Pw EFFECT ON THE TANG PART OF CONTACT

      }
      KRATOS_CATCH( "" )

   }



   // TANGENT MATRIX OF THE TANGENT CONTACT. It has to be better programed (but it works).
   void AxisymPointRigidContactPenaltyWater2DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      MatrixType InitialMatrix = rLeftHandSideMatrix;
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;
      Matrix  WaterMatrix = ZeroMatrix(MatSize,MatSize);

      // CHECK THE YIELD CRITERION
      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

      double TangentForceModulus = this->CalculateCoulombsFrictionLaw( rVariables.Gap.Tangent, NormalForceModulus, rVariables); // ALSO COMPUTES the TANGENT MATRIX ( and is saved in the rVariables.TangentMatrix.Normal and .Tangent)
      // OBS: rVariables.TangentMatrix.Normal is the variation of the tangent stress with respect to the normal stress
      //      rVariables.TangentMatrix.Normal is the variation of the tangent stress with respect the tangent gap

      if( fabs(TangentForceModulus) >= 1e-25 ){

         // PLASTIC
         if ( rVariables.Slip ) {


            if ( rVariables.Penalty.Tangent < 1e-6)
               return ;

            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
            Matrix  TangentMatrix = ZeroMatrix(dimension,dimension);

            TangentMatrix = rVariables.TangentMatrix.Normal * ( rVariables.Penalty.Normal / rVariables.ContributoryFactor) * outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Normal) ;
            TangentMatrix += rVariables.TangentMatrix.Tangent * outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Tangent);

            TangentMatrix *= rIntegrationWeight * rVariables.ContributoryFactor;
            for (unsigned int i = 0; i < MatSize-1; i++){
               for (unsigned int j = 0; j < MatSize-1; j++){
                  rLeftHandSideMatrix(i,j) -= TangentMatrix(i,j);
               }
            }

            // and now the water pressure part
            for (unsigned int i = 0; i < MatSize-1; i++) {
               WaterMatrix(i, MatSize-1) += rVariables.TangentMatrix.Normal * rVariables.Surface.Tangent(i) * rIntegrationWeight * rVariables.ContributoryFactor ;

            }
            rLeftHandSideMatrix -= WaterMatrix;
         }
         // ELASTIC
         else {
            Matrix SmallMatrix =  rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent) ;
            for (unsigned int i = 0; i < MatSize-1; i++){
               for (unsigned int j = 0; j < MatSize-1; j++){
                  rLeftHandSideMatrix(i,j) += SmallMatrix(i,j);
               }
            }


         }
      }

      KRATOS_CATCH( "" )

      return;



   }
   //***********************************************************************************
   //***********************************************************************************
   // CALCULATE CONTACT FORCE, the same but then adds a zero on the water row.
   void AxisymPointRigidContactPenaltyWater2DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)

   {
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

      VectorType  LocalForce = ZeroVector(MatSize-1);
      AxisymPointRigidContactPenalty2DCondition::CalculateAndAddContactForces( LocalForce, rVariables, rIntegrationWeight);

      for (unsigned int i = 0; i < MatSize-1; i++){
         rRightHandSideVector(i) = LocalForce(i);
      }
      rRightHandSideVector(MatSize-1) = 0;


   }


} // end namespace kratos
