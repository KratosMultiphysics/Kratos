//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//


#include "custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_wP_3D_condition.hpp"


// vale, això és allò però amb el terme de les pressions d'aigua
// he de copiar bé els dofs dels elements i tal

namespace Kratos
{
   // USUAL CONSTRUCTORS

   //************************************************************************************
   //************************************************************************************
   EPPointRigidContactPenaltywP3DCondition::EPPointRigidContactPenaltywP3DCondition(IndexType NewId, GeometryType::Pointer
         pGeometry)
      : EPPointRigidContactPenalty3DCondition(NewId, pGeometry)
   {
   }

   //************************************************************************************
   //************************************************************************************
   EPPointRigidContactPenaltywP3DCondition::EPPointRigidContactPenaltywP3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : EPPointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties)
   {
   }


   //************************************************************************************
   //************************************************************************************
   EPPointRigidContactPenaltywP3DCondition::EPPointRigidContactPenaltywP3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : EPPointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties, pRigidWall)
   {

   }

   //************************************************************************************
   //************************************************************************************
   EPPointRigidContactPenaltywP3DCondition::EPPointRigidContactPenaltywP3DCondition( EPPointRigidContactPenaltywP3DCondition const& rOther )
      : EPPointRigidContactPenalty3DCondition(rOther)
   {
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer EPPointRigidContactPenaltywP3DCondition::Create(IndexType NewId, NodesArrayType
         const& ThisNodes,  PropertiesType::Pointer pProperties) const
   {
     return Kratos::make_shared<EPPointRigidContactPenaltywP3DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
   }

   Condition::Pointer EPPointRigidContactPenaltywP3DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
   {
      EPPointRigidContactPenaltywP3DCondition NewCondition( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
      NewCondition.mCurrentInfo = this->mCurrentInfo;
      NewCondition.mSavedInfo   = this->mSavedInfo;

      return Kratos::make_shared<EPPointRigidContactPenaltywP3DCondition>(NewCondition);

   }


   //************************************************************************************
   //************************************************************************************
   EPPointRigidContactPenaltywP3DCondition::~EPPointRigidContactPenaltywP3DCondition()
   {
   }

   // AQUI VIENEN LAS COSAS

   // ADD A NEW DOF DUE TO THE WATER
   void EPPointRigidContactPenaltywP3DCondition::GetDofList(DofsVectorType& rConditionDofList,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      rConditionDofList.resize(0);
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
         if ( dimension == 3)
            rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
         rConditionDofList.push_back(GetGeometry()[i].pGetDof(WATER_PRESSURE));
      }


      KRATOS_CATCH( "" )
   }

   //***********************************************************************************
   //***********************************************************************************

   void EPPointRigidContactPenaltywP3DCondition::EquationIdVector(EquationIdVectorType& rResult,
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
         if (dimension == 3) {
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
         } else {
            rResult[index + 2] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
         }


      }

      KRATOS_CATCH( "" )
   }


   //***********************************************************************************
   //***********************************************************************************

   void EPPointRigidContactPenaltywP3DCondition::GetValuesVector(Vector& rValues, int Step)
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
         if (dimension == 3) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
         } else {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
         }
      }
   }

   //***********************************************************************************
   //***********************************************************************************

   void EPPointRigidContactPenaltywP3DCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
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
         if ( dimension == 3) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = 0.0;
         } else {
            rValues[index + 2] = 0.0;
         }

      }
   }


   //***********************************************************************************
   //***********************************************************************************

   void EPPointRigidContactPenaltywP3DCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
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
         if (dimension == 3) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[index + 3] = 0.0;
         } else {
            rValues[index + 2] = 0.0;
         }
      }

   }


   // MODIFY THE INIZIALIZATION due to the differnt size of the system due to the water pressure
   void EPPointRigidContactPenaltywP3DCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * ( dimension + 1) ; // water pressure

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

   //****************** CalculateAndAddKuug (since matrices are of different size) ************
   //******************************************************************************************

   void EPPointRigidContactPenaltywP3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      if( rVariables.Options.Is(ACTIVE)){

         Matrix KuugN(3,3);
         noalias(KuugN) = ZeroMatrix(3,3);

         noalias(KuugN) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);

         Matrix KuugT(4,4);
         noalias(KuugT) = ZeroMatrix(4,4);

         this->CalculateAndAddKuugTangent( KuugT,  rVariables, rIntegrationWeight);


         for(unsigned int i=0; i<dimension; i++)
         {
            for(unsigned int j=0; j<dimension; j++)
            {
               rLeftHandSideMatrix(i,j) += KuugN(i,j) + KuugT(i,j);
            }
            rLeftHandSideMatrix(i, dimension) += KuugT(i,3);
         }


      }
      else{
         rLeftHandSideMatrix= ZeroMatrix(dimension+1,dimension+1);
      }

      // if( rVariables.Options.Is(ACTIVE))
      //   std::cout<<" Contact Tangent Matrix ["<<this->Id()<<"]: "<<rLeftHandSideMatrix<<std::endl;

      KRATOS_CATCH( "" )
   }

   // ********************* CalculateKuugTangent (with effective water pressure term) *******************
   // ***************************************************************************************************
   void EPPointRigidContactPenaltywP3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

      ConstitutiveVariables ConstVariables;
      Vector AuxVector;

      rVariables.Slip = this->CalculateFrictionLaw( rVariables, ConstVariables, AuxVector);
      double TangentForceModulus = norm_2( AuxVector);

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix  LHSMatrix = ZeroMatrix(dimension, dimension);
      double Area = this->CalculateSomeSortOfArea();


      if( fabs(TangentForceModulus) >= 1e-25 ){

        MatrixType Identity(3,3);
        noalias(Identity) = IdentityMatrix(3);

        if( rVariables.Slip ){

          noalias(LHSMatrix) += ConstVariables.TangentTangentMatrix * outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) * rIntegrationWeight * Area;
          noalias(LHSMatrix) -= ConstVariables.NormalTangentMatrix * outer_prod( ConstVariables.ForceDirection, rVariables.Surface.Normal) * rVariables.Penalty.Normal * rIntegrationWeight;
          noalias(LHSMatrix) += ConstVariables.TangentForceRatio * rVariables.Penalty.Tangent * ( Identity - outer_prod( rVariables.Surface.Normal, rVariables.Surface.Normal) - outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) ) * rIntegrationWeight ;

        }
        else {


          noalias(LHSMatrix) = rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent); // aquest terme jo no el veig

          noalias(LHSMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * ( Identity - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

        }

      }


      // ensamble all togheter
      for (unsigned int i = 0; i < dimension; i++) {
        for (unsigned int j = 0; j < dimension; j++) {
          rLeftHandSideMatrix(i,j) += LHSMatrix(i,j);
        }
      }

      if ( rVariables.Slip && ( fabs(TangentForceModulus) >= 1e-25)) {
        Vector LHSVector = ZeroVector(3);
        noalias( LHSVector ) = ConstVariables.NormalTangentMatrix * ConstVariables.ForceDirection * rIntegrationWeight * Area ;
        for (unsigned int i = 0; i < dimension; i++) {
          rLeftHandSideMatrix(i, dimension) += LHSVector(i);
        }
      }

      KRATOS_CATCH("")
   }


} // end namespace kratos
