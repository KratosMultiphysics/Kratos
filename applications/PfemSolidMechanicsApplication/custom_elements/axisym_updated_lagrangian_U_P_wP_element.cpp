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
#include "custom_elements/axisym_updated_lagrangian_U_P_wP_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"
//#include "includes/define.h"
//#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"


// U P wP formulation with added water pressure

// there are some geometric-like terms related to the water pressure missing
//
// the darcy law geometric terms are missing


namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   AxisymUpdatedLagrangianUPwPElement::AxisymUpdatedLagrangianUPwPElement()
      : AxisymUpdatedLagrangianUPressureElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPwPElement::AxisymUpdatedLagrangianUPwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : AxisymUpdatedLagrangianUPressureElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPwPElement::AxisymUpdatedLagrangianUPwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : AxisymUpdatedLagrangianUPressureElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPwPElement::AxisymUpdatedLagrangianUPwPElement( AxisymUpdatedLagrangianUPwPElement const& rOther)
      :AxisymUpdatedLagrangianUPressureElement(rOther)
       , mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPwPElement&  AxisymUpdatedLagrangianUPwPElement::operator=(AxisymUpdatedLagrangianUPwPElement const& rOther)
   {
      AxisymUpdatedLagrangianUPressureElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUPwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUPwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUPwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUPwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new AxisymUpdatedLagrangianUPwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPwPElement::~AxisymUpdatedLagrangianUPwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void AxisymUpdatedLagrangianUPwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );


         rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
      }
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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

         rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
         rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::GetValuesVector( Vector& rValues, int Step )
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

         rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
         rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
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
         rValues[index + 2] = 0;
         rValues[index + 3] = 0;
      }
   }

   //*********************************ACCELERATION***************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
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

         rValues[index + 2] = 0;
         rValues[index + 3] = 0;
      }

   }


   void AxisymUpdatedLagrangianUPwPElement::InitializeElementData ( ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      AxisymUpdatedLagrangianUPressureElement::InitializeElementData( rVariables, rCurrentProcessInfo );

      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

   }

   //************************************************************************************
   //************************************************************************************

   int  AxisymUpdatedLagrangianUPwPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = AxisymUpdatedLagrangianUPressureElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P Second element type ", " AxisymUpdatedLagrangianUPwPElement" )

            //verify that the variables are correctly initialized

            if ( WATER_PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

      if ( this->GetProperties().Has(THICKNESS) ) {
	      double thickness = this->GetProperties()[THICKNESS];
	      if ( thickness <= 0.0) {
		      this->GetProperties()[THICKNESS] = 1.0;
	      }
      } else {
	     this->GetProperties()[THICKNESS] = 1.0;
      } 
                  return correct;

      KRATOS_CATCH( "" );
   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
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

   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();


      // I HAVE TO DO SOMETHING WITH THE GENERAL VARIABLES * DETFT ETC
      LocalSystemComponents UPLocalSystem; 
      unsigned int MatSize = number_of_nodes * ( dimension+1);
      MatrixType  LocalLeftHandSideMatrix = ZeroMatrix(MatSize,MatSize);

      UPLocalSystem.SetLeftHandSideMatrix( LocalLeftHandSideMatrix);


      AxisymUpdatedLagrangianUPressureElement::CalculateAndAddLHS( UPLocalSystem, rVariables, rIntegrationWeight);


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
      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      CalculateAndAddUnconsideredKuuTerms( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKuwP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPu( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPwP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPwPStab( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;


   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      /*if ( this->Id() == 1) {
        std::cout <<  " ID1 rHS RHS rhs R: detF0 " << rVariables.detF0 << " detF " << rVariables.detF << std::endl;
        }*/

      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
      CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, IntegrationWeight );

      // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
      CalculateAndAddInternalForces( rRightHandSideVector, rVariables, IntegrationWeight);

      // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
      CalculateAndAddPressureForces( rRightHandSideVector, rVariables, IntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, IntegrationWeight);

      // operation performed rRight
      CalculateAndAddWaterPressureForces( rRightHandSideVector, rVariables, IntegrationWeight );

      CalculateAndAddStabilizedWaterPressure( rRightHandSideVector, rVariables, IntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
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

   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int voigtsize = 4;

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

      TotalStress(0) = StressVector(0);
      TotalStress(1) = StressVector(1);
      TotalStress(2) = StressVector(2);
      TotalStress(3) = StressVector(3);



      double WaterPressure = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];

      for (unsigned int i = 0; i < 3; i++) {
         TotalStress(i) += WaterPressure; 
      }

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), TotalStress );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexup = ( dimension + 2 ) * i ;
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
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddWaterPressureForces( VectorType& rRightHandSideVector, 
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      VectorType Fh=rRightHandSideVector;


      Matrix K = ZeroMatrix(dimension,dimension);
      Vector b = ZeroVector(dimension);
      double WaterDensity = GetProperties()[DENSITY_WATER];
      b(dimension-1) = -10.0*WaterDensity;
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      unsigned int indexp = dimension + 1;
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

            rRightHandSideVector[indexp] += (1.0/WaterBulk) * rVariables.N[i] * rVariables.N[j] * DeltaPressure * rIntegrationWeight * ScalingConstant / rVariables.detF0;

            for ( unsigned int p = 0; p < dimension; ++p )
            {

               rRightHandSideVector[indexp] -= rVariables.N[i]*DeltaDisplacement[p] * rVariables.DN_DX(j,p) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               if (p == 0)
                  rRightHandSideVector[indexp] -= rVariables.N[i]*DeltaDisplacement[p] * rVariables.N[j]*(1.0/ rVariables.CurrentRadius) * rIntegrationWeight * ScalingConstant/ rVariables.detF0;

               for ( unsigned int q = 0; q < dimension; ++q )
               {

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

      KRATOS_CATCH( "" )

   }

   // ******************************* CALCULATE AND ADD PRESSURE FORCES ***************
   // *********************************************************************************
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      VectorType Fh=rRightHandSideVector;

      double ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3; i++)
         ElementalMeanStress += rVariables.StressVector(i) / 3.0;


      unsigned int indexp = dimension;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {


            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber; //2D

         }

         rRightHandSideVector[indexp] -= ElementalMeanStress * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;


         indexp += (dimension + 2);
      }


      KRATOS_CATCH( "" )

   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


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
         LargeDisplacementElement::GetValueOnIntegrationPoints( SHEAR_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization /= Values[0];

         LargeDisplacementElement::GetValueOnIntegrationPoints( BULK_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization *= Values[0];

      }


      unsigned int indexp = dimension;
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

         }


         indexp += (dimension + 2);
      }



      KRATOS_CATCH( "" )

   }


   // ********* STABILIZATION OF THE MASS BALANCE EQUATION **************************
   // *******************************************************************************
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddStabilizedWaterPressure( VectorType& rRightHandSideVector, 
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant ; 
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      double StabilizationAlpha, Caux, StabilizationFactor; 

      const ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      LargeDisplacementElement::GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      if ( Mmodulus[0] < 0.0001) {
         LargeDisplacementElement::GetValueOnIntegrationPoints( YOUNG_MODULUS, Mmodulus, CurrentProcessInfo); 
         Mmodulus[0] = GetProperties()[YOUNG_MODULUS];
      }
      Caux = 1.0/Mmodulus[0];

      double he;
      he = GetElementSize( rVariables.DN_DX);
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];


      StabilizationAlpha = pow(he, 2) * Caux / (6.0) - DeltaTime * Permeability / 2.0;

      if (StabilizationAlpha < 0.0)
      {
         return;
      }

      //StabilizationAlpha *= 0.5  + 0.5 * tanh( pow(he, 2) * Caux / (Permeability* DeltaTime)  );
      StabilizationAlpha *= StabilizationFactor;

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      VectorType Fh=rRightHandSideVector;


      Matrix K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = 1.0;

      unsigned int indexp = dimension+1;

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
      // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )




   }

   // ** ***************** Calculation of the geometric terms due to the water pressure 
   //
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddUnconsideredKuuTerms(MatrixType& rK,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // ************** Calculation of the Ku wP Matrix *********************************************
   //
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddKuwP(MatrixType& rLeftHandSideMatrix,
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
         unsigned int indexup = (dimension +2) * i ;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;

               if (k==0)
                  rLeftHandSideMatrix(indexup+k, indexp) +=rVariables.N(i) * rVariables.N[j] * (1.0/ rVariables.CurrentRadius)* rIntegrationWeight* rVariables.detF;
            }
            indexp += (dimension + 2);
         }
      }

      // std::cout<<std::endl;
      // std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }

   //* *************************** Calculation of the KwP U Matrix
   //
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddKwPu(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      MatrixType Kh=rLeftHandSideMatrix;

      unsigned int indexp = dimension+1;

      for ( unsigned int i = 0; i < number_of_nodes; ++i )
      {
         for ( unsigned int j = 0; j < number_of_nodes; ++j )
         {
            int indexup = (dimension+2)*j ;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexp, indexup+k) += rVariables.N[i] * rVariables.DN_DX( j , k ) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               if (k==0)
                  rLeftHandSideMatrix(indexp, indexup+k) +=rVariables.N[i]*rVariables.N[j] * rIntegrationWeight*ScalingConstant* (1.0/ rVariables.CurrentRadius) / rVariables.detF0;

               // TRY TO PROGRAM LD TERMS...
               /*for (unsigned int m = 0; m < number_of_nodes; ++m ) 
               {
                  const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[m].FastGetSolutionStepValue( DISPLACEMENT );
                  const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[m].FastGetSolutionStepValue( DISPLACEMENT , 1 );
                  array_1d<double, 3 > DeltaDisplacement            = CurrentDisplacement-PreviousDisplacement;
                  const double & WaterPressure = GetGeometry()[m].FastGetSolutionStepValue( WATER_PRESSURE );
                  for (unsigned n = 0; n < dimension; n++) {

                     rLeftHandSideMatrix(indexp, indexup+k) -= rVariables.N[i] * DeltaDisplacement[n] * rVariables.DN_DX(m,k) * rVariables.DN_DX(j,n) * rIntegrationWeight * ScalingConstant / rVariables.detF0;

                     for (unsigned t = 0; t < dimension; ++t) {
                        rLeftHandSideMatrix(indexp, indexup+k) += (rVariables.DN_DX(i,k)*rVariables.DN_DX(m,n) + rVariables.DN_DX(i,n)*rVariables.DN_DX(m,k) ) * K(t,n) * rVariables.DN_DX(j,t) * WaterPressure * DeltaTime * rIntegrationWeight * ScalingConstant / rVariables.detF0;
                     }
                  }
               } */

            }
         }
         indexp += (dimension + 2);
      }



      // std::cout<<std::endl;
      // std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )

   }

   // ********************** Calculation of the KwP P Matrix
   // 
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddKwPP(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      // this is zero!!!!

      KRATOS_CATCH("")
   }

   //** ***************** Calculation of the K wP wP Matrix
   //
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddKwPwP(MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight )
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);
      Matrix K = ZeroMatrix(dimension,dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension+1;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension+1;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(WaterBulk)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * ScalingConstant / rVariables.detF0;
            for ( unsigned int p = 0; p < dimension; ++p ) {
               for ( unsigned int q = 0; q < dimension; ++q ) {
                  rLeftHandSideMatrix(indexpi, indexpj) -= DeltaTime * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q) * rIntegrationWeight* ScalingConstant / rVariables.detF0;
               }
            }

            indexpj += (dimension + 2);
         }

         indexpi += (dimension + 2);
      }

      // std::cout<<std::endl;
      // std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) <<std::endl;
      // std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) / ScalingConstant / rIntegrationWeight* rVariables.detF0  * WaterBulk <<std::endl;
      KRATOS_CATCH( "" )

   }

   // ** ***************** Calculation of the Stabilization Tangent Matrix
   //
   void AxisymUpdatedLagrangianUPwPElement::CalculateAndAddKwPwPStab(MatrixType& rLeftHandSideMatrix,
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

      const ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      LargeDisplacementElement::GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      if ( Mmodulus[0] < 0.0001) {
         //GetValueOnIntegrationPoints( YOUNG_MODULUS, Mmodulus, CurrentProcessInfo); 
         Mmodulus[0] = GetProperties()[YOUNG_MODULUS];
      }
      Caux = 1.0/Mmodulus[0];

      double he;
      he = GetElementSize( rVariables.DN_DX);
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];


      // nova Manera
      StabilizationAlpha = pow(he, 2) * Caux / (6.0) - DeltaTime * Permeability / 2.0;

      if ( StabilizationAlpha < 0.0)
         return;

      // StabilizationAlpha *= 0.5  + 0.5 * tanh( pow(he, 2) * Caux / (Permeability* DeltaTime)  );
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
      unsigned int indexpi = dimension+1;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension+1;
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
      //std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) / ScalingConstant / rIntegrationWeight* rVariables.detF0  * WaterBulk <<std::endl;


      KRATOS_CATCH( "" )


   }
   //************************************CALCULATE VOLUME CHANGE*************************
   //************************************************************************************

   double& AxisymUpdatedLagrangianUPwPElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
   {
      KRATOS_TRY

      rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

      return rVolumeChange;

      KRATOS_CATCH( "" )
   }

   ////************************************************************************************
   ////************************************************************************************


   // GET THE (GEOMETRICAL) SIZE FOR THE STABILIZATION TERM
   // this size is recomended from Sun, Ostien and Salinger (IJNAMG, 2013)

   double AxisymUpdatedLagrangianUPwPElement::GetElementSize( const Matrix& rDN_DX)
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
   void AxisymUpdatedLagrangianUPwPElement::GetConstants(double& rScalingConstant, double& rWaterBulk, double& rDeltaTime, double& rPermeability)
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

   void AxisymUpdatedLagrangianUPwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUPressureElement )
   }

   void AxisymUpdatedLagrangianUPwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUPressureElement )
   }


}
