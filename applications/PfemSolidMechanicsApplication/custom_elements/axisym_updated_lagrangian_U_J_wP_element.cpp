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

#include "custom_utilities/axisym_water_pressure_utilities_Jacobian.hpp"

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


   //************* GETTING METHODS ******************************************************
   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {

      KRATOS_TRY


      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rOutput.size() != integration_points_number )
         rOutput.resize( integration_points_number );
      /*if ( rVariable == EFF_CON_WATER_FORCE) {

         // vale, torno a fer de les meves...
         ElementDataType Variables;
         this->InitializeElementData( Variables, rCurrentProcessInfo);

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
         ElementDataType Variables;
         this->InitializeElementData( Variables, rCurrentProcessInfo);

         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions=Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
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
            this->SetElementData(Variables,Values,PointNumber);

            // OBS, now changing Variables I change Values because they are pointers ( I hope);
            double ElementalDetFT = Variables.detH;
            Matrix ElementalFT = Variables.H;

            // AND NOW IN THE OTHER WAY
            Matrix m; double d; 
            ComputeConstitutiveVariables( Variables, m, d);

            Variables.H = m;
            Variables.detH = d; 
            Values.SetDeformationGradientF( Variables.H);
            Values.SetDeterminantF( Variables.detH );

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            Variables.H = ElementalFT;
            Variables.detH = ElementalDetFT;


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
      else {*/
         UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      //}

      KRATOS_CATCH("")
   }



   void AxisymUpdatedLagrangianUJwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
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
         rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
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
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
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
            rValues[index + 2] = 0;
            rValues[index + 3] = 0;
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
            rValues[index + 2] = 0;
            rValues[index + 3] = 0;
         }

   }


   // ********************** Initialize General Variables ***********************************
   // ***************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::InitializeElementData ( ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      AxisymUpdatedLagrangianUJElement::InitializeElementData( rVariables, rCurrentProcessInfo );

      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

      //stabilization factor
     double StabilizationFactor = 1.0;
     if( GetProperties().Has(STABILIZATION_FACTOR_WP) ){
       StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];
     }
     else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR_WP) ){
       StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR_WP];
     }
     GetProperties().SetValue(STABILIZATION_FACTOR_WP, StabilizationFactor);

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
         KRATOS_THROW_ERROR( std::invalid_argument, "WATER PRESSURE has Key zero! (check if the application is correctly registered", "" )

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


   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      KRATOS_TRY

      // define some variables
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // ComputeBaseClass LHS
      LocalSystemComponents UJLocalSystem; 
      unsigned int MatSize = number_of_nodes * ( dimension+1);
      MatrixType  LocalLeftHandSideMatrix = ZeroMatrix(MatSize,MatSize) ;
      UJLocalSystem.SetLeftHandSideMatrix( LocalLeftHandSideMatrix);

      // LHS. base class
      AxisymUpdatedLagrangianUJElement::CalculateAndAddLHS( UJLocalSystem, rVariables, rIntegrationWeight);

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      // Reshape the BaseClass LHS and Add the Hydro Part
      AxisymWaterPressureJacobianUtilities WaterUtility;
      Matrix TotalF = prod ( rVariables.F, rVariables.F0);
      int number_of_variables = dimension + 2; // displ - Jacobian - waterPressure
      Vector VolumeForce;
      VolumeForce= this->CalculateVolumeForce( VolumeForce, rVariables);
// 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetDeformationGradient( TotalF);
      HMVariables.SetVolumeForce( VolumeForce);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      HMVariables.CurrentRadius = rVariables.CurrentRadius;
      //HMVariables.ConstrainedModulus
      HMVariables.number_of_variables = number_of_variables;


      // LHS. hydromechanical problem
     rLeftHandSideMatrix = WaterUtility.CalculateAndAddHydromechanicalLHS( HMVariables, rLeftHandSideMatrix, LocalLeftHandSideMatrix, IntegrationWeight);

      // LHS. water pressure stabilization
      ProcessInfo SomeProcessInfo; 
      std::vector< double> Mmodulus;
      GetValueOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }
HMVariables.ConstrainedModulus = ConstrainedModulus;
      rLeftHandSideMatrix = WaterUtility.CalculateAndAddStabilizationLHS( HMVariables, rLeftHandSideMatrix, IntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      KRATOS_TRY

      // define some variables
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;


      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // Compute Base Class RHS
      LocalSystemComponents BaseClassLocalSystem;
      Vector BaseClassRightHandSideVector = ZeroVector ( (dimension+1) * number_of_nodes );
      BaseClassLocalSystem.SetRightHandSideVector(BaseClassRightHandSideVector );
      Vector VolumeForce = rVolumeForce; 
      VolumeForce *= 0.0;
      AxisymUpdatedLagrangianUJElement::CalculateAndAddRHS( BaseClassLocalSystem, rVariables, VolumeForce, rIntegrationWeight);

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      // Reshape the BaseClass RHS and Add the Hydro Part
      AxisymWaterPressureJacobianUtilities WaterUtility;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);

      int number_of_variables = dimension+2; // displacement - Jacobian - waterPressure

      // 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetDeformationGradient( TotalF);
      HMVariables.SetVolumeForce( rVolumeForce);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      HMVariables.CurrentRadius = rVariables.CurrentRadius;
      //HMVariables.ConstrainedModulus
      HMVariables.number_of_variables = number_of_variables;

      rRightHandSideVector = WaterUtility.CalculateAndAddHydromechanicalRHS( HMVariables, rRightHandSideVector, BaseClassRightHandSideVector, IntegrationWeight); 

      // Add Stab term
      ProcessInfo SomeProcessInfo; 
      std::vector<double> Mmodulus;
      GetValueOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }
      HMVariables.ConstrainedModulus = ConstrainedModulus;
      rRightHandSideVector = WaterUtility.CalculateAndAddStabilization( HMVariables, rRightHandSideVector, IntegrationWeight); 

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;


      KRATOS_CATCH( "" )


   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUJElement )
   }

   void AxisymUpdatedLagrangianUJwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUJElement )
   }


}
