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
#include "custom_elements/axisym_updated_lagrangian_U_J_wP_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"
//#include "includes/define.h"
//#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"


namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   AxisymUpdatedLagrangianUJwPElement::AxisymUpdatedLagrangianUJwPElement()
      : AxisymUpdatedLagrangianUJElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJwPElement::AxisymUpdatedLagrangianUJwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : AxisymUpdatedLagrangianUJElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJwPElement::AxisymUpdatedLagrangianUJwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : AxisymUpdatedLagrangianUJElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJwPElement::AxisymUpdatedLagrangianUJwPElement( AxisymUpdatedLagrangianUJwPElement const& rOther)
      :AxisymUpdatedLagrangianUJElement(rOther)
       , mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJwPElement&  AxisymUpdatedLagrangianUJwPElement::operator=(AxisymUpdatedLagrangianUJwPElement const& rOther)
   {
      AxisymUpdatedLagrangianUJElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUJwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUJwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new AxisymUpdatedLagrangianUJwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJwPElement::~AxisymUpdatedLagrangianUJwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************
   void AxisymUpdatedLagrangianUJwPElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
   }
   void AxisymUpdatedLagrangianUJwPElement::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
   }
   void AxisymUpdatedLagrangianUJwPElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {

      KRATOS_TRY


      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rOutput.size() != integration_points_number )
         rOutput.resize( integration_points_number );
      if ( rVariable == EFF_CON_WATER_FORCE) {

         // vale, torno a fer de les meves...
         GeneralVariables Variables;
         this->InitializeGeneralVariables( Variables, rCurrentProcessInfo);

         const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod);
         //reading integration points
         for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            const unsigned int number_of_nodes = GetGeometry().PointsNumber();

            Vector StressVector = ZeroVector(4); // we are axisym

            double WaterPressure = 0;
            for (unsigned int i = 0; i < number_of_nodes; i++)
               WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * Variables.N[i];

            for (unsigned int i = 0; i < 3; i++) {
               StressVector(i) += WaterPressure; 
            }

            double thisInt = integration_points[PointNumber].Weight() * Variables.detJ;
            double IntegrationWeight = thisInt * 2.0 * 3.151692654 * Variables.CurrentRadius; 
            Vector WPF = IntegrationWeight * prod( trans( Variables.B ), StressVector );

            rOutput[PointNumber] = WPF; 
         }

      }
      else if ( ( rVariable == EFF_CON_TOTAL_FORCE) || ( rVariable == EFF_CON_EFFEC_FORCE) )
      {
         // vale, torno a fer de les meves...
         GeneralVariables Variables;
         this->InitializeGeneralVariables( Variables, rCurrentProcessInfo);

         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions=Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

         const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod);
         //reading integration points

         
         for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {

            // 1. COmpute The Cauchy stress

            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //to take in account previous step writing
            if( mFinalizedStep ){
               this->GetHistoricalVariables(Variables,PointNumber);
            }		

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            // OBS, now changing Variables I change Values because they are pointers ( I hope);
            double ElementalDetFT = Variables.detFT;
            Matrix ElementalFT = Variables.FT;

            // AND NOW IN THE OTHER WAY
            Matrix m; double d; 
            ComputeConstitutiveVariables( Variables, m, d);

            Variables.FT = m;
            Variables.detFT = d; 
            Values.SetDeformationGradientF( Variables.FT);
            Values.SetDeterminantF( Variables.detFT );

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            Variables.FT = ElementalFT;
            Variables.detFT = ElementalDetFT;


            // 2. Do My THIngs
            Vector StressVector = Variables.StressVector; 

            const unsigned int number_of_nodes = GetGeometry().PointsNumber();
            // 2.b Compute The total stress (if necessary) 
            if ( rVariable == EFF_CON_TOTAL_FORCE )
            {

               double WaterPressure = 0.0;
               for (unsigned int i = 0; i < number_of_nodes; i++)
                  WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * Variables.N[i];

               for (unsigned int i = 0; i < 3; i++) {
                  StressVector(i) += WaterPressure; 
               }
            }

            double thisInt = integration_points[PointNumber].Weight() * Variables.detJ;
            double IntegrationWeight = thisInt * 2.0 * 3.151692654 * Variables.CurrentRadius; 
            Vector WPF = IntegrationWeight * prod( trans( Variables.B ), StressVector );

            rOutput[PointNumber] = WPF; 
         }


      }
      else {
         UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }

      KRATOS_CATCH("")
   }



   void AxisymUpdatedLagrangianUJwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

         rElementalDofList.push_back( GetGeometry()[i].pGetDof( JACOBIAN ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
      }
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
            rResult[index + 3] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }
         else
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::GetValuesVector( Vector& rValues, int Step )
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
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }
         else
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
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

   void AxisymUpdatedLagrangianUJwPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
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


   void AxisymUpdatedLagrangianUJwPElement::InitializeGeneralVariables ( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      AxisymUpdatedLagrangianUJElement::InitializeGeneralVariables( rVariables, rCurrentProcessInfo );

      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

   }

   //************************************************************************************
   //************************************************************************************

   int  AxisymUpdatedLagrangianUJwPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = AxisymUpdatedLagrangianUJElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-J element type ", " AxisymUpdatedLagrangianUJwPElement" )

            //verify that the variables are correctly initialized

            if ( WATER_PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }




   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
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

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();


      // I HAVE TO DO SOMETHING WITH THE GENERAL VARIABLES * DETFT ETC
      LocalSystemComponents UJLocalSystem; 
      unsigned int MatSize = number_of_nodes * ( dimension+1);
      MatrixType  LocalLeftHandSideMatrix = ZeroMatrix(MatSize,MatSize) ;

      UJLocalSystem.SetLeftHandSideMatrix( LocalLeftHandSideMatrix);


      AxisymUpdatedLagrangianUJElement::CalculateAndAddLHS( UJLocalSystem, rVariables, rIntegrationWeight);


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

      CalculateAndAddUnconsideredKuuTerms( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      CalculateAndAddKuwP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPu( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPJ( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPwP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      CalculateAndAddKwPwPStab( rLeftHandSideMatrix, rVariables, IntegrationWeight);

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

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
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
      CalculateAndAddJacobianForces( rRightHandSideVector, rVariables, IntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedJacobian( rRightHandSideVector, rVariables, IntegrationWeight);

      // operation performed rRight
      CalculateAndAddWaterPressureForces( rRightHandSideVector, rVariables, IntegrationWeight );

      CalculateAndAddStabilizedWaterPressure( rRightHandSideVector, rVariables, IntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
         GeneralVariables& rVariables,
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

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      VectorType Fh=rRightHandSideVector;

      Vector StressVector = rVariables.StressVector;

      double WaterPressure = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];

      for (unsigned int i = 0; i < 3; i++) {
         StressVector(i) += WaterPressure; 
      }

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );

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

   //******************************** JACOBIAN FORCES  **********************************
   //************************************************************************************
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            const double& NodalJac = (GetGeometry()[j].GetSolutionStepValue( JACOBIAN) );

            rRightHandSideVector[indexp] +=   rVariables.N[i] * rVariables.N[j]  * NodalJac * rIntegrationWeight / rVariables.detFT;

         }

         rRightHandSideVector[indexp] -= rVariables.N[i] * rVariables.detFT * rIntegrationWeight / rVariables.detFT;



         indexp += (dimension + 2);

      }

      KRATOS_CATCH( "" )

   }



   //****************** STABILIZATION *********************************************************
   //************************* defined in the Stab element ************************************

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
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
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
      AlphaStabilization *= StabilizationFactor; 

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

      AlphaStabilization *= BulkModulus;  // TIMES THE BULK MODULUS BECAUSE I HAVE ALL THE EQUATION MULTIPLIED BY THE BULK MODULUS

      if (YoungModulus < 0.00001)
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0;
         AlphaStabilization *= mElementStabilizationNumber; 

      }

      double consistent = 1;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            consistent=(-1.0)*AlphaStabilization;
            if(i==j)
               consistent=2.0*AlphaStabilization;

            double& Jacobian= GetGeometry()[j].FastGetSolutionStepValue(JACOBIAN);
            rRightHandSideVector[indexp] += consistent * Jacobian * rIntegrationWeight / (rVariables.detFT);

            // std::cout<<" Pressure "<<Pressure<<std::endl;
         }


         indexp += (dimension + 2);
      }


      // std::cout<<std::endl;
      // std::cout<<" IntegrationWeight "<<rIntegrationWeight<<" detF "<<rVariables.detF0<<std::endl;
      // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )

   }

   // ********* MASS BALANCE EQUATION: WATER PRESSURE EQUATION ***********************
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddWaterPressureForces( VectorType& rRightHandSideVector, 
         GeneralVariables& rVariables,
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


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {


            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            const double& CurrentNodalJacobian = GetGeometry()[j].FastGetSolutionStepValue( JACOBIAN );
            const double& PreviousNodalJacobian = GetGeometry()[j].FastGetSolutionStepValue( JACOBIAN , 1 );
            double DeltaNodalJacobian = CurrentNodalJacobian - PreviousNodalJacobian; 


            // water compressibility term
            rRightHandSideVector[indexp] += (1.0/WaterBulk) * rVariables.N[i] * rVariables.N[j] * DeltaPressure * rIntegrationWeight * ScalingConstant / rVariables.detF0;

            // mixture volume change.
            rRightHandSideVector[indexp] -= rVariables.N[i] * rVariables.N[j] * DeltaNodalJacobian / ElementalNodalJacobian * rIntegrationWeight * ScalingConstant / rVariables.detF0;

            for ( unsigned int p = 0; p < dimension; ++p )
            {

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

   // ********* STABILIZATION OF THE MASS BALANCE EQUATION **************************
   // *******************************************************************************
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddStabilizedWaterPressure( VectorType& rRightHandSideVector, 
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant ; 
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      double StabilizationAlpha, Caux, StabilizationFactor; 

      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      if ( Mmodulus[0] < 0.0001) {
         GetValueOnIntegrationPoints( YOUNG_MODULUS, Mmodulus, CurrentProcessInfo); 
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


      unsigned int indexp = dimension + 1;
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
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddUnconsideredKuuTerms(MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      return; 
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      int size = number_of_nodes * dimension;

      Matrix Kuu = zero_matrix<double>(size,size);

      // axisymmetric geometric matrix

      double alpha1 = 0;
      double alpha2 = 0;
      double alpha3 = 0;


      Vector StressVector = ZeroVector(4);
      double WaterPressure = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];

      for (unsigned int i = 0;  i < 3; i++)
         StressVector(i) = WaterPressure; 


      unsigned int indexi = 0;
      unsigned int indexj = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         indexj =0;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            alpha1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * StressVector(0) + rVariables.DN_DX(i,1) * StressVector(3) );
            alpha2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * StressVector(3) + rVariables.DN_DX(i,1) * StressVector(1) );
            alpha3 = rVariables.N[i] * rVariables.N[j] * StressVector(2) * (1.0/rVariables.CurrentRadius*rVariables.CurrentRadius);

            Kuu(indexi,indexj)     = alpha1 + alpha2 + alpha3 ;
            Kuu(indexi+1,indexj+1) = alpha1 + alpha2 ;

            indexj+=2;
         }

         indexi+=2;

      }

      Kuu *= rIntegrationWeight;

      //std::cout<<std::endl;
      //std::cout<<" Kuu "<<Kuu<<std::endl;


      MatrixType Kh=rLeftHandSideMatrix;

      //assemble into rk the geometric uu contribution:
      indexi = 0;
      indexj = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int idim = 0; idim < dimension ; idim ++)
         {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
               {
                  rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                  indexj++;
               }
            }
            indexi++;
         }
      }

      //std::cout<<std::endl;
      //std::cout<<" Kgeo "<<rK-Kh<<std::endl;

      KRATOS_CATCH( "" )

   }

   // ************** Calculation of the Ku wP Matrix *********************************************
   //
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddKuwP(MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
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
               if ( k == 0)
                  rLeftHandSideMatrix(indexup+k, indexp) += rVariables.N[i] * rVariables.N[j] * ( 1.0/ rVariables.CurrentRadius) * rIntegrationWeight * rVariables.detF; 
            }
            indexp += (dimension + 2);
         }
      }

      KRATOS_CATCH("")
   }

   //* *************************** Calculation of the KwP U Matrix
   //
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddKwPu(MatrixType& rK,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // ********************** Calculation of the KwP J Matrix
   // 
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddKwPJ(MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      Matrix SmallMatrix = ZeroMatrix(number_of_nodes,number_of_nodes);

      double ElementalNodalJacobian = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         ElementalNodalJacobian += GetGeometry()[i].FastGetSolutionStepValue( JACOBIAN) * rVariables.N[i];

      double consistent; 
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            consistent = 1.0/12.0;
            if ( i == j)
               consistent *= 2.0;

            SmallMatrix(i,j) += consistent / ElementalNodalJacobian;
            for (unsigned int k = 0; k < number_of_nodes; k++) {

               consistent = 1.0/12.0;
               if ( i == k)
                  consistent *= 2.0;
               const double& CurrentNodalJacobian = GetGeometry()[k].FastGetSolutionStepValue( JACOBIAN );
               const double& PreviousNodalJacobian = GetGeometry()[k].FastGetSolutionStepValue( JACOBIAN , 1 );
               double DeltaJacobian = CurrentNodalJacobian- PreviousNodalJacobian; 

               SmallMatrix(i,j) -= consistent * DeltaJacobian / pow(ElementalNodalJacobian, 2) * rVariables.N[j] ;
            }
         }
      }




      SmallMatrix *= rIntegrationWeight * ScalingConstant / rVariables.detF0;

      // AND NOW I HAVE TO PUT THE MATRIX IN ITS POSITION

      unsigned int dime = dimension + 2;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            rLeftHandSideMatrix(dime*(i+1)-1, dime*(j+1)-2) += SmallMatrix(i,j); // .
         }
      }



      KRATOS_CATCH("")
   }

   //** ***************** Calculation of the K wP wP Matrix
   //
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddKwPwP(MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
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

      Matrix SmallMatrix = ZeroMatrix(number_of_nodes,number_of_nodes);

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

      for (unsigned int i = 0; i< number_of_nodes; i++) {
         for (unsigned int j = 0; j< number_of_nodes; j++) {
            SmallMatrix(i,j) += ( 1.0/WaterBulk) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * ScalingConstant/ rVariables.detF0;
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
   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddKwPwPStab(MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
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
      GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
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
         unsigned int indexpj = dimension +1;
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

   double& AxisymUpdatedLagrangianUJwPElement::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
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

   double AxisymUpdatedLagrangianUJwPElement::GetElementSize( const Matrix& rDN_DX)
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
   void AxisymUpdatedLagrangianUJwPElement::GetConstants(double& rScalingConstant, double& rWaterBulk, double& rDeltaTime, double& rPermeability)
   {
      //double ScalingConstant = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
      rScalingConstant = 1.0;
      double Something = GetProperties()[MY_SCALING_CONSTANT];
      if ( Something > 1.0e-6) {
         rScalingConstant = Something; 
      }
      rDeltaTime  = GetProperties()[DELTA_TIME];
      rDeltaTime = mTimeStep;
      rPermeability = GetProperties()[PERMEABILITY];
      rWaterBulk = GetProperties()[WATER_BULK_MODULUS];

   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void AxisymUpdatedLagrangianUJwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }







}
