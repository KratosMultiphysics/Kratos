//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                 -0.1 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_P_wP_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"

// U wP formulation with added water pressure

// there are some geometric-like terms related to the water pressure missing
//
// the darcy law geometric terms are missing

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUPwPElement::UpdatedLagrangianUPwPElement()
      : UpdatedLagrangianUPressureElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPwPElement::UpdatedLagrangianUPwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUPressureElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPwPElement::UpdatedLagrangianUPwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUPressureElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUPwPElement::UpdatedLagrangianUPwPElement( UpdatedLagrangianUPwPElement const& rOther)
      :UpdatedLagrangianUPressureElement(rOther)
       , mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUPwPElement&  UpdatedLagrangianUPwPElement::operator=(UpdatedLagrangianUPwPElement const& rOther)
   {
      UpdatedLagrangianUPressureElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUPwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUPwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUPwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUPwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

      //-----------//

      if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
         NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

      for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
      {
         NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
      }

      NewElement.mDeterminantF0 = mDeterminantF0;

      NewElement.mTimeStep = mTimeStep; 

      return Element::Pointer( new UpdatedLagrangianUPwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPwPElement::~UpdatedLagrangianUPwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void UpdatedLagrangianUPwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

         rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
      }
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rResult.size() != element_size )
         rResult.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int index = i * ( dimension + 2);
         rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

         if( dimension == 3)
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }
         else
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::GetValuesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension +2) ;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }
         else
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension + 2);  
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = 0;
            rValues[index + 4] = 0;
         }
         else
         {
            rValues[index + 2] = 0;
            rValues[index + 3] = 0;
         }
      }
   }

   //*********************************ACCELERATION***************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension + 2);
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = 0;
            rValues[index + 4] = 0;
         }
         else
         {
            rValues[index + 2] = 0;
            rValues[index + 3] = 0;
         }
      }

   }


   void UpdatedLagrangianUPwPElement::InitializeElementData ( ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      UpdatedLagrangianUPressureElement::InitializeElementData( rVariables, rCurrentProcessInfo );

      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

   }

   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUPwPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = UpdatedLagrangianUPressureElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P Second element type ", " UpdatedLagrangianUPwPElement" )

            //verify that the variables are correctly initialized

            if ( WATER_PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rCalculationFlags.Is(LargeDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
         if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

         noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }


      //resizing as needed the RHS
      if ( rCalculationFlags.Is(LargeDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
         if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

         rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
      }
   }

   //************************************************************************************
   //************************************************************************************


   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************




   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();


      // I HAVE TO DO SOMETHING WITH THE GENERAL VARIABLES * DETFT ETC
      LocalSystemComponents UPLocalSystem; 
      unsigned int MatSize = number_of_nodes * ( dimension+1);
      MatrixType  LocalLeftHandSideMatrix = ZeroMatrix(MatSize, MatSize) ;

      UPLocalSystem.SetLeftHandSideMatrix( LocalLeftHandSideMatrix);


      UpdatedLagrangianUPressureElement::CalculateAndAddLHS( UPLocalSystem, rVariables, rIntegrationWeight);


      // PUT THE MATRIX IT IN THE RIGHT POSITION ( the waterPressure dofs are in the middle of the system matrix )
      unsigned int dime1 = dimension+1;
      unsigned int dime2 = dimension+2;
      for (unsigned int iNode = 0; iNode < number_of_nodes; iNode++ ) {
         for (unsigned int idim = 0; idim < dimension +1; idim++ ) {
            for (unsigned int jNode = 0; jNode < number_of_nodes; jNode++ ) {
               for (unsigned int jdim = 0; jdim < dimension+1; jdim++) {
                  rLeftHandSideMatrix( dime2*iNode + idim, dime2*jNode+jdim) += LocalLeftHandSideMatrix( dime1*iNode + idim, dime1*jNode + jdim);
               }
            } 
         }
      }

      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      CalculateAndAddUnconsideredKuuTerms( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      CalculateAndAddKuwP( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      CalculateAndAddKwPu( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      CalculateAndAddKwPP( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      CalculateAndAddKwPwP( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      CalculateAndAddKwPwPStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      if ( this->Id() ==0)
      {
         std::cout << " SMALL MATRIX " << LocalLeftHandSideMatrix << std::endl;
         std::cout << " TANGENT MATRIX " << rLeftHandSideMatrix << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
      }

   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      /*if ( this->Id() == 1) {
        std::cout <<  " ID1 rHS RHS rhs R: detF0 " << rVariables.detF0 << " detF " << rVariables.detF << std::endl;
        }*/

      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
      CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

      // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
      CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
      CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

      // operation performed rRight
      CalculateAndAddWaterPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight );

      CalculateAndAddStabilizedWaterPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
         ElementDataType& rVariables,
         Vector& rVolumeForce,
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //VectorType Fh=rRightHandSideVector;

      double DomainChange = (1.0/rVariables.detF0); //density_n+1 = density_0 * ( 1.0 / detF0 )

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int indexup = ( dimension + 2) * i;
         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] += rIntegrationWeight * rVariables.N[i] * rVolumeForce[j] * DomainChange;
         }

      }

      // std::cout<<std::endl;
      // std::cout<<" Fext "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )
   }


   //************************** INTERNAL FORCES    *******************************
   //************************************** Idem but with Total Stress ***********

   void UpdatedLagrangianUPwPElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int voigtsize = 3;
      if (dimension == 3) 
         voigtsize = 6;

      VectorType Fh=rRightHandSideVector;

      Vector StressVector = rVariables.StressVector;

      double ElemMeanPressure = 0.0; 
      for (unsigned int i = 0; i < 3; i++)
         ElemMeanPressure += StressVector(i);
      ElemMeanPressure /= 3.0;
      double NodalMeanPressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         NodalMeanPressure += GetGeometry()[i].GetSolutionStepValue( PRESSURE ) * rVariables.N[i];

      for (unsigned int i = 0; i < 3; i++)
         StressVector(i) += ( NodalMeanPressure - ElemMeanPressure );

      Vector TotalStress = ZeroVector(voigtsize); 

      if ( voigtsize == 6) 
      {
         TotalStress = StressVector; 
      }
      else {
         TotalStress(0) = StressVector(0);
         TotalStress(1) = StressVector(1);
         TotalStress(2) = StressVector(3);
      }



      double WaterPressure = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];

      for (unsigned int i = 0; i < dimension; i++) {
         TotalStress(i) += WaterPressure; 
      }

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), TotalStress );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexup = ( dimension +2 ) * i ;
         unsigned int indexu  = dimension * i;

         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
         }
      }

      KRATOS_CATCH( "" )
   }



   // ********* MASS BALANCE EQUATION: WATER PRESSURE EQUATION ***********************
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   void UpdatedLagrangianUPwPElement::CalculateAndAddWaterPressureForces( VectorType& rRightHandSideVector, 
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension + 1;

      VectorType Fh=rRightHandSideVector;


      Vector b = ZeroVector(dimension);
      double WaterDensity = GetProperties()[DENSITY_WATER];
      b(dimension-1) = -10.0*WaterDensity;

      Matrix K;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);
      //this->GetPermeabilityTensor( Permeability, TotalF, K);
      K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; i++)
         K(i,i) = Permeability;

      double ElementalNodalJacobian = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         ElementalNodalJacobian += GetGeometry()[i].FastGetSolutionStepValue( JACOBIAN) * rVariables.N[i];


      double consistent;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[j].FastGetSolutionStepValue( DISPLACEMENT );
            const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[j].FastGetSolutionStepValue( DISPLACEMENT , 1 );
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;


            consistent = 1.0/12.0;
            if ( i == j)
               consistent = 2.0/12.0;

            // water compressibility term
            rRightHandSideVector[indexp] += (1.0/WaterBulk) * consistent * DeltaPressure * rIntegrationWeight * ScalingConstant / rVariables.detF0;


            for ( unsigned int p = 0; p < dimension; ++p )
            {
               // mixture volume change.
               rRightHandSideVector[indexp] -= rVariables.N[i]*DeltaDisplacement[p] * rVariables.DN_DX(j,p) * rIntegrationWeight * ScalingConstant / rVariables.detF0;

               for ( unsigned int q = 0; q < dimension; ++q )
               {
                  // DARCY FLOW
                  rRightHandSideVector[indexp] += DeltaTime * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q)* CurrentPressure *rIntegrationWeight * ScalingConstant / rVariables.detF0;

                  if (j == 0) {
                     // Gravity term of the DARCY FLOW.
                     rRightHandSideVector[indexp] += DeltaTime * rVariables.DN_DX(i,p) * K(p,q)*b(q)*rIntegrationWeight* ScalingConstant / rVariables.detF0;
                  }

               }


            }

         }


         indexp += (dimension + 2);

      }

      // std::cout<<std::endl;
      // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
      // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH("")
   }

   // ******************************* CALCULATE AND ADD PRESSURE FORCES ***************
   // *********************************************************************************
   void UpdatedLagrangianUPwPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;
      VectorType Fh=rRightHandSideVector;
      double consistent=1;

      double ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3; i++)
         ElementalMeanStress += rVariables.StressVector(i) / 3.0;


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            consistent=1;
            if(i==j)
               consistent=2;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += consistent * (1.0/12.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber; //2D

         }

         rRightHandSideVector[indexp] -= ElementalMeanStress * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;


         indexp += (dimension + 2);
      }


      KRATOS_CATCH( "" )

   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      // VectorType Fh=rRightHandSideVector;
      // std::cout<<" Element "<<this->Id()<<" "<<std::endl;

      //use of this variable for the complete parameter: (deffault: 4)
      double AlphaStabilization  = 4.0; 
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_P];
      AlphaStabilization *= StabilizationFactor;

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));
      AlphaStabilization *= BulkModulus;


      if ( YoungModulus < 0.00001 )
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0 ;

         ProcessInfo SomeProcessInfo;
         std::vector<double> Values;

         SolidElement::GetValueOnIntegrationPoints( SHEAR_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization /= Values[0];

         SolidElement::GetValueOnIntegrationPoints( BULK_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization *= Values[0];

      }


      double consistent = 1;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            consistent=(-1)*AlphaStabilization;
            if(i==j)
               consistent=2*AlphaStabilization;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;

            // std::cout<<" Pressure "<<Pressure<<std::endl;
         }


         indexp += (dimension + 2);
      }



      KRATOS_CATCH( "" )

   }


   // ********* STABILIZATION OF THE MASS BALANCE EQUATION **************************
   // *******************************************************************************
   void UpdatedLagrangianUPwPElement::CalculateAndAddStabilizedWaterPressure( VectorType& rRightHandSideVector, 
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {
      // FPL just copied from there
      KRATOS_TRY

      double ScalingConstant ;
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      double StabilizationAlpha, Caux, StabilizationFactor; 

      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;

      SolidElement::GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      Caux = 1.0/Mmodulus[0];

      double he;
      he = GetElementSize( rVariables.DN_DX);
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];

      StabilizationAlpha = pow(he, 2.0) * Caux / ( 6.0) - DeltaTime*Permeability; 
      StabilizationAlpha *= StabilizationFactor;

      if (StabilizationAlpha < 0.0)
      {
         return;
      }

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      Matrix K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = 1.0;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {


            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            for ( unsigned int p = 0; p < dimension; ++p )
            {

               for ( unsigned int q = 0; q < dimension; ++q )
               {

                  rRightHandSideVector[indexp] += StabilizationAlpha * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q)* DeltaPressure *rIntegrationWeight * ScalingConstant / rVariables.detF0;

               }


            }

         }

         indexp += (dimension + 2);

      }


      // std::cout<<std::endl;
      // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
      // std::cout<<" Fpres STABILI "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )


   }

   // ** ***************** Calculation of the geometric terms due to the water pressure 
   //
   void UpdatedLagrangianUPwPElement::CalculateAndAddUnconsideredKuuTerms(MatrixType& rK,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // ************** Calculation of the Ku wP Matrix *********************************************
   //
   void UpdatedLagrangianUPwPElement::CalculateAndAddKuwP(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //MatrixType Kh=rLeftHandSideMatrix;
      //contributions to stiffness matrix calculated on the reference configuration
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexp  = dimension+1;
         unsigned int indexup = (dimension +2 )  * i ;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;
            }
            indexp += (dimension + 2);
         }
      }

      KRATOS_CATCH("")
   }

   //* *************************** Calculation of the KwP U Matrix
   //
   void UpdatedLagrangianUPwPElement::CalculateAndAddKwPu(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix K;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);
      //this->GetPermeabilityTensor( Permeability, TotalF, K);
      K = ZeroMatrix(3,3);
      for (unsigned int i = 0; i < 3; i++)
         K(i,i) = Permeability;

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexp = dimension+1;


      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            int indexup = (dimension+2)*j ;
            for (unsigned int k = 0; k < dimension; k++)
            {
               // Small Strains solid skeleton volume change
               rLeftHandSideMatrix(indexp, indexup + k) += rVariables.N[i] * rVariables.DN_DX( j, k) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
            }
         }
         indexp += (dimension + 2);
      }


      // std::cout<<std::endl;
      // std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH("")
   }

   // ********************** Calculation of the KwP P Matrix
   // 
   void UpdatedLagrangianUPwPElement::CalculateAndAddKwPP(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY


      KRATOS_CATCH("")
   }

   //** ***************** Calculation of the K wP wP Matrix
   //
   void UpdatedLagrangianUPwPElement::CalculateAndAddKwPwP(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight )
   {
      KRATOS_TRY

      // mmmm
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      Matrix K;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);
      //this->GetPermeabilityTensor( Permeability, TotalF, K);
      K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; i++)
         K(i,i) = Permeability;

      Matrix SmallMatrix = ZeroMatrix( number_of_nodes, number_of_nodes);

      for (unsigned int m = 0; m < number_of_nodes; m++) {   
         for (unsigned int n = 0; n < number_of_nodes; n++) {   
            for (unsigned int i = 0; i < dimension; i++) {   
               for (unsigned int j = 0; j < dimension; j++) {  
                  SmallMatrix(m,n) +=  rVariables.DN_DX(m,i) * K(i,j) * rVariables.DN_DX(n,j) ;
               }
            }
         }
      }
      SmallMatrix *= DeltaTime * rIntegrationWeight * ScalingConstant / rVariables.detF0;

      double consistent ;
      for (unsigned int i = 0; i< number_of_nodes; i++) {
         for (unsigned int j = 0; j< number_of_nodes; j++) {
            consistent = 1.0/12.0;
            if ( i == j) 
               consistent *= 2.0;
            SmallMatrix(i,j) += ( 1.0/WaterBulk) * consistent * rIntegrationWeight * ScalingConstant/ rVariables.detF0;
         }
      }

      unsigned int dime = dimension + 2;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            rLeftHandSideMatrix(dime*(i+1)-1, dime*(j+1)-1) -= SmallMatrix(i,j);
         }
      }



      KRATOS_CATCH("")
   }

   // ** ***************** Calculation of the Stabilization Tangent Matrix
   //
   void UpdatedLagrangianUPwPElement::CalculateAndAddKwPwPStab(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      double StabilizationAlpha, Caux, StabilizationFactor; 

      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      SolidElement::GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      Caux = 1.0/Mmodulus[0];

      double he;
      he = GetElementSize( rVariables.DN_DX);
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];

      StabilizationAlpha = pow(he, 2.0) * Caux / ( 6.0) - DeltaTime*Permeability; 
      StabilizationAlpha *= StabilizationFactor;

      if (StabilizationAlpha < 0.0)
      {
         return;
      }

      Matrix K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = 1.0;

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            for ( unsigned int p = 0; p < dimension; ++p ) {
               for ( unsigned int q = 0; q < dimension; ++q ) {
                  rLeftHandSideMatrix(indexpi, indexpj) -= StabilizationAlpha * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q) * rIntegrationWeight* ScalingConstant / rVariables.detF0;
               }
            }

            indexpj += (dimension + 2);
         }

         indexpi += (dimension + 2);
      }

      // std::cout<<std::endl;
      // std::cout<<" Kpp STABILI "<< (rLeftHandSideMatrix-Kh) <<std::endl;
      //std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) / ScalingConstant / rIntegrationWeight* rVariables.detF0  * WaterBulk <<std::endl;


      KRATOS_CATCH("")
   }
   //************************************CALCULATE VOLUME CHANGE*************************
   //************************************************************************************

   double& UpdatedLagrangianUPwPElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
   {
      KRATOS_TRY

      rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

      return rVolumeChange;

      KRATOS_CATCH( "" )
   }

   ////************************************************************************************
   ////************************************************************************************

   void UpdatedLagrangianUPwPElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      //create and initialize element variables:
      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      //create constitutive law parameters:
      ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

      //set constitutive law flags:
      Flags &ConstitutiveLawOptions=Values.GetOptions();

      ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
      ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


      for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
      {

         //compute element kinematics B, F, DN_DX ...
         this->CalculateKinematics(Variables,PointNumber);

         //set general variables to constitutivelaw parameters
         this->SetElementData(Variables,Values,PointNumber);


         //call the constitutive law to update material variables
         mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse(Values, Variables.StressMeasure);

         //call the constitutive law to finalize the solution step
         mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
               GetGeometry(),
               Variables.N,
               rCurrentProcessInfo );



         //call the element internal variables update
         this->FinalizeStepVariables(Variables,PointNumber);
      }


      this->Set(SolidElement::FINALIZED_STEP,true);

      KRATOS_CATCH( "" )
   }
   ////************************************************************************************
   ////************************************************************************************

   ////************************************************************************************
   ////************************************************************************************

   ////************************************************************************************
   ////************************************************************************************

   void UpdatedLagrangianUPwPElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      //create and initialize element variables:
      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      //create constitutive law parameters:
      ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

      //set constitutive law flags:
      Flags &ConstitutiveLawOptions=Values.GetOptions();

      ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
      ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
      ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

      //reading integration points
      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      //auxiliary terms
      Vector VolumeForce;

      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
         //compute element kinematics B, F, DN_DX ...
         this->CalculateKinematics(Variables,PointNumber);

         //set general variables to constitutivelaw parameters
         this->SetElementData(Variables,Values,PointNumber);


         //compute stresses and constitutive parameters
         mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, Variables.StressMeasure);


         //some transformation of the configuration can be needed (UL element specially)
         this->TransformElementData(Variables,PointNumber);

         //calculating weights for integration on the "reference configuration"
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
         IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


         if ( rLocalSystem.CalculationFlags.Is(LargeDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
         {
            //contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );
         }

         if ( rLocalSystem.CalculationFlags.Is(LargeDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
         {
            //contribution to external forces
            VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );

            this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, IntegrationWeight );
         }
      }


      KRATOS_CATCH( "" )
   }

   // GET THE (GEOMETRICAL) SIZE FOR THE STABILIZATION TERM
   // this size is recomended from Sun, Ostien and Salinger (IJNAMG, 2013)

   double UpdatedLagrangianUPwPElement::GetElementSize( const Matrix& rDN_DX)
   {
      double he = 0.0;

      unsigned int number_of_nodes = rDN_DX.size1();
      unsigned int dimension = rDN_DX.size2();

      double aux;
      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         aux = 0;
         for (unsigned int p = 0; p < dimension ; p++)
         {
            aux += rDN_DX(i,p);
         }
         he += fabs(aux);
      }
      he *= sqrt( double(dimension) );
      he = 4.0/he;

      return he;

   }


   // SOMETHING STUPID
   void UpdatedLagrangianUPwPElement::GetConstants(double& rScalingConstant, double& rWaterBulk, double& rDeltaTime, double& rPermeability)
   {
      //double ScalingConstant = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
      rScalingConstant = 1.0;
      rDeltaTime  = GetProperties()[DELTA_TIME];
      rDeltaTime = mTimeStep;
      rPermeability = GetProperties()[PERMEABILITY];
      rWaterBulk = GetProperties()[WATER_BULK_MODULUS];

   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUPressureElement )
   }

   void UpdatedLagrangianUPwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUPressureElement )
   }


}
