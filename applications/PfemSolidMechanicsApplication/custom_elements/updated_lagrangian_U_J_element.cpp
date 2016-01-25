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
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_U_J_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

#include "pfem_solid_mechanics_application.h"

// MIRAR POR ENCIMA DE IMPLEMENTAR UNA FORMULACION U-J, de esta forma cuando llamo a la constitutiva pongo Fbarra de los desplazamientos y el jacobiano de J.
// Como a la parte constitutiva solo paso punteros (i no matrices y otras cosas) se vuelve un poco asín
//
// THE RESIDUALS ARE OK
// TO DO: Improve the tangent matrices (I think now they are correct).
// VALE, PER MI QUE TOT ESTÀ MÉS O MENYS BÉ.

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUJElement::UpdatedLagrangianUJElement()
      : LargeDisplacementElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJElement::UpdatedLagrangianUJElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : LargeDisplacementElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJElement::UpdatedLagrangianUJElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : LargeDisplacementElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUJElement::UpdatedLagrangianUJElement( UpdatedLagrangianUJElement const& rOther)
      :LargeDisplacementElement(rOther)
       ,mDeformationGradientF0(rOther.mDeformationGradientF0)
       ,mDeterminantF0(rOther.mDeterminantF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUJElement&  UpdatedLagrangianUJElement::operator=(UpdatedLagrangianUJElement const& rOther)
   {
      LargeDisplacementElement::operator=(rOther);

      mDeformationGradientF0.clear();
      mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
      }

      mDeterminantF0 = rOther.mDeterminantF0;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUJElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUJElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUJElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJElement::~UpdatedLagrangianUJElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void UpdatedLagrangianUJElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

         rElementalDofList.push_back( GetGeometry()[i].pGetDof( JACOBIAN ));
      }
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = number_of_nodes * dimension + number_of_nodes;

      if ( rResult.size() != element_size )
         rResult.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int index = i * dimension + i;
         rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

         if( dimension == 3)
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
         }
         else
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
         }

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::GetValuesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension + i;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
         }
         else
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension + i;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = 0;
         }
         else
         {
            rValues[index + 2] = 0;
         }
      }
   }

   //*********************************ACCELERATION***************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension + i;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = 0;
         }
         else
         {
            rValues[index + 2] = 0;
         }
      }

   }

   // **************************************************************************
   // **************************************************************************
   void UpdatedLagrangianUJElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      LargeDisplacementElement::InitializeSolutionStep( rCurrentProcessInfo );

      // VALE, AQUÍ POSO LES COSES ESTUPIDES QUE JO VULGUI
      const double& YoungModulus = GetProperties()[YOUNG_MODULUS];
      if (YoungModulus < 0.00001)
      {
         std::vector<double> Values;
         GetValueOnIntegrationPoints( SIMILAR_SHEAR_MODULUS, Values, rCurrentProcessInfo);
         mElementStabilizationNumber = 1.0 / Values[0];

         GetValueOnIntegrationPoints( SIMILAR_BULK_MODULUS, Values, rCurrentProcessInfo);
         mElementStabilizationNumber *= Values[0];


      }

      KRATOS_CATCH( "" );
   }

   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUJElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-J element type ", " UpdatedLagrangianUJElement" )

            //verify that the variables are correctly initialized

            if ( PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }

   //*********************************SET DOUBLE VALUE***********************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
         std::vector<double>& rValues,
         const ProcessInfo& rCurrentProcessInfo )
   {

      if (rVariable == DETERMINANT_F){

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();


         for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
         {
            mDeterminantF0[PointNumber] = rValues[PointNumber];
         }

      }
      else{

         LargeDisplacementElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }


   }


   //**********************************GET DOUBLE VALUE**********************************
   //************************************************************************************


   void UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
         std::vector<double>& rValues,
         const ProcessInfo& rCurrentProcessInfo )
   {

      const unsigned int& integration_points_number = mConstitutiveLawVector.size();
      if (rVariable == DETERMINANT_F){


         if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

         for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
         {
            rValues[PointNumber] = mDeterminantF0[PointNumber];
         }

      }
      else{

         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }

   }

   void UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {

      KRATOS_TRY


      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rOutput.size() != integration_points_number )
         rOutput.resize( integration_points_number );
      if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
      {
         //create and initialize element variables:
         GeneralVariables Variables;
         this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions=Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

         //reading integration points
         for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {
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



            //call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR)
               mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
            else
               mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);

            Variables.FT = ElementalFT;
            Variables.detFT = ElementalDetFT;

            if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
               rOutput[PointNumber].resize( Variables.StressVector.size(), false );

            rOutput[PointNumber] = Variables.StressVector;


         }

      }
      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }

      KRATOS_CATCH("")
   }


   void UpdatedLagrangianUJElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == CONSTITUTIVE_MATRIX )
      {
         //create and initialize element variables:
         GeneralVariables Variables;
         this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions=Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

         //reading integration points
         for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

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

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            Variables.FT = ElementalFT;
            Variables.detFT = ElementalDetFT;


            if( rOutput[PointNumber].size2() != Variables.ConstitutiveMatrix.size2() )
               rOutput[PointNumber].resize( Variables.ConstitutiveMatrix.size1() , Variables.ConstitutiveMatrix.size2() , false );

            rOutput[PointNumber] = Variables.ConstitutiveMatrix;

         }


      }
      else if (rVariable == EQ_CAUCHY_STRESS) 
      {
         //create and initialize element variables:
         GeneralVariables Variables;
         this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

         Variables.StressVector = ZeroVector(6); // I WANT TO GET THE THIRD COMPONENT
         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions=Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

         //reading integration points
         for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {
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

            if ( ( rOutput[PointNumber].size1() != 3 ) |
                  ( rOutput[PointNumber].size2() != 3 ) )
               rOutput[PointNumber].resize( 3, 3, false );

            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor( Variables.StressVector );


         }

      }
      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }
   }

   void UpdatedLagrangianUJElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {


      if ( rVariable == DETERMINANT_F) {
         const unsigned int& integration_points_number = mConstitutiveLawVector.size();

         if (rOutput.size() != integration_points_number)
            rOutput.resize( integration_points_number) ;

         for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
         {
            rOutput[PointNumber] = mDeterminantF0[PointNumber];
         }
      }
      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }

   }



   void UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( const Variable<Vector> & rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
   {

      LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
   }

   //**********************************GET TENSOR VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == EQ_CAUCHY_STRESS) {
         CalculateOnIntegrationPoints(rVariable, rValue, rCurrentProcessInfo);
      }
      else if ( rVariable == CAUCHY_STRESS_VECTOR ) {
         CalculateOnIntegrationPoints( EQ_CAUCHY_STRESS, rValue, rCurrentProcessInfo);
      }
      else {
         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);
      }

   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************
   void UpdatedLagrangianUJElement::Initialize()
   {
      KRATOS_TRY

      LargeDisplacementElement::Initialize();

      SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //Resize historic deformation gradient
      if ( mDeformationGradientF0.size() != integration_points_number )
         mDeformationGradientF0.resize( integration_points_number );

      if ( mDeterminantF0.size() != integration_points_number )
         mDeterminantF0.resize( integration_points_number, false );

      for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
         mDeterminantF0[PointNumber] = 1;
         mDeformationGradientF0[PointNumber] = identity_matrix<double> (dimension);
      }

      const unsigned int number_of_nodes = GetGeometry().size();
      for (unsigned int Node = 0; Node < number_of_nodes; Node++) {
         double& DetFNodal = GetGeometry()[Node].GetSolutionStepValue(JACOBIAN );
         DetFNodal = 1.0;
      }

      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      LargeDisplacementElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

      //Calculate Delta Position
      rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

      //set variables including all integration points values

      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );


   }

   ////************************************************************************************
   ////************************************************************************************

   void UpdatedLagrangianUJElement::FinalizeStepVariables( GeneralVariables & rVariables, const double& rPointNumber )
   { 
      //update internal (historical) variables
      mDeterminantF0[rPointNumber]         = rVariables.detF * rVariables.detF0 ;
      mDeformationGradientF0[rPointNumber] = prod( rVariables.F, rVariables.F0 ) ;
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

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


   void UpdatedLagrangianUJElement::CalculateDeformationMatrix(Matrix& rB,
         Matrix& rF,
         Matrix& rDN_DX)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      rB.clear(); //set all components to zero

      if( dimension == 2 )
      {

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

         }

      }
      else if( dimension == 3 )
      {

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );

         }
      }
      else
      {

         KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

      }

      KRATOS_CATCH( "" )
   }



   //*************************COMPUTE DEFORMATION GRADIENT*******************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateDeformationGradient(const Matrix& rDN_DX,
         Matrix& rF,
         Matrix& rDeltaPosition)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      rF = identity_matrix<double> ( dimension );

      if( dimension == 2 )
      {

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
         }

      }
      else if( dimension == 3)
      {

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {

            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 0 , 2 ) += rDeltaPosition(i,0)*rDN_DX ( i , 2 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
            rF ( 1 , 2 ) += rDeltaPosition(i,1)*rDN_DX ( i , 2 );
            rF ( 2 , 0 ) += rDeltaPosition(i,2)*rDN_DX ( i , 0 );
            rF ( 2 , 1 ) += rDeltaPosition(i,2)*rDN_DX ( i , 1 );
            rF ( 2 , 2 ) += rDeltaPosition(i,2)*rDN_DX ( i , 2 );
         }

      }
      else
      {

         KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

      }

      KRATOS_CATCH( "" )
   }

   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************


   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************
   void UpdatedLagrangianUJElement::CalculateKinematics(GeneralVariables& rVariables,
         const double& rPointNumber)

   {
      KRATOS_TRY

      //Get the parent coodinates derivative [dN/d£]
      const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

      //Get the shape functions for the order of the integration method [N]
      const Matrix& Ncontainer = rVariables.GetShapeFunctions();

      //Parent to reference configuration
      rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
      Matrix InvJ;
      MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

      //Compute cartesian derivatives [dN/dx_n]
      noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

      //Current Deformation Gradient [dx_n+1/dx_n]
      //this->CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition);

      //Deformation Gradient F [dx_n+1/dx_n] to be updated
      noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

      //Determinant of the deformation gradient F
      rVariables.detF  = MathUtils<double>::Det(rVariables.F);

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
      Matrix Invj;
      MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

      //Compute cartesian derivatives [dN/dx_n+1]
      rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

      //Determinant of the Deformation Gradient F0
      rVariables.detF0 = mDeterminantF0[rPointNumber];
      rVariables.F0    = mDeformationGradientF0[rPointNumber];

      //Set Shape Functions Values for this integration point
      rVariables.N=row( Ncontainer, rPointNumber);

      //Compute the deformation matrix B
      this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
   {


      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // operation performed: add Km to the rLefsHandSideMatrix

      //respect to the current configuration n+1
      CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kg to the rLefsHandSideMatrix
      CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kup to the rLefsHandSideMatrix
      CalculateAndAddKuJ( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpu to the rLefsHandSideMatrix
      CalculateAndAddKJu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpp to the rLefsHandSideMatrix
      CalculateAndAddKJJ( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpp Stab to the rLefsHandSideMatrix
      CalculateAndAddKJJStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;


   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      /*if ( this->Id() == 1) {
        std::cout <<  " ID1 rHS RHS rhs R: detF0 " << rVariables.detF0 << " detF " << rVariables.detF << std::endl;
        }*/
      // EM FALTA UN CACHO
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
      CalculateAndAddJacobianForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedJacobian( rRightHandSideVector, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
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
         int indexup = dimension * i + i;
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

   void UpdatedLagrangianUJElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      VectorType Fh=rRightHandSideVector;

      Vector StressVector = rVariables.StressVector;

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexup = dimension * i + i;
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
   void UpdatedLagrangianUJElement::CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      double consistent; 
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent = 1.0/12.0;
            if ( i == j)
               consistent *= 2.0;

            const double& NodalJac = (GetGeometry()[j].GetSolutionStepValue( JACOBIAN) );

            rRightHandSideVector[indexp] +=   consistent  * NodalJac * rIntegrationWeight / rVariables.detFT;

         }

         rRightHandSideVector[indexp] -= rVariables.N[i] * rVariables.detFT * rIntegrationWeight / rVariables.detFT;



         indexp += (dimension + 1);

      }

      KRATOS_CATCH( "" )

   }



   //****************** STABILIZATION *********************************************************
   //************************* defined in the Stab element ************************************

   void UpdatedLagrangianUJElement::CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
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
            rRightHandSideVector[indexp] += consistent * Jacobian * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

            // std::cout<<" Pressure "<<Pressure<<std::endl;
         }


         indexp += (dimension + 1);
      }


      // std::cout<<std::endl;
      // std::cout<<" IntegrationWeight "<<rIntegrationWeight<<" detF "<<rVariables.detF0<<std::endl;
      // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )

   }


   //******** Kuu Material************************************************************
   //***************** It includes the pw geometric stiffness ************************

   void UpdatedLagrangianUJElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      double dimension_double = double(dimension);

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;

      unsigned int voigtsize = 3;
      if (dimension == 3) 
         voigtsize = 6;

      Matrix DeviatoricTensor = ZeroMatrix(voigtsize);
      Vector Identity = ZeroVector(voigtsize);
      for (unsigned int i = 0; i < voigtsize ; ++i) {
         DeviatoricTensor(i,i) = 1.0;
      }
      for (unsigned int i = 0; i < dimension; i++) {
         Identity(i) = 1.0;
         for (unsigned int j = 0; j < dimension; j++) {
            DeviatoricTensor( i,j) -= 1.0/dimension_double;
         }
      }


      ConstitutiveMatrix = prod( ConstitutiveMatrix, DeviatoricTensor);

      if ( this->Id() == 0) 
      {
         std::cout << " CONS 0 " << rVariables.ConstitutiveMatrix << std::endl;
         std::cout << " CONS 1 " << ConstitutiveMatrix << std::endl;
      }
      Matrix AuxMatrix = ZeroMatrix(voigtsize);

      for (unsigned int i = 0; i < voigtsize; i++) {
         for (unsigned int j = 0; j < voigtsize; j++) {
            ConstitutiveMatrix(i,j)  += (1 -  2/dimension_double) * rVariables.StressVector(i) * Identity(j);
            AuxMatrix(i,j) += rVariables.StressVector(i) * Identity(j);
         }
      }

      if ( this->Id() == 0) 
      {
         std::cout << " CONS 2 " << ConstitutiveMatrix << std::endl;
         std::cout << " AUX MATRIX " << AuxMatrix << std::endl;
         std::cout << " stressVector " << rVariables.StressVector << std::endl;
         std::cout << std::endl;
      }


      Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( ConstitutiveMatrix, rVariables.B ) ) ); 

      // MatrixType Kh=rLeftHandSideMatrix;

      unsigned int indexi = 0;
      unsigned int indexj  = 0;
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

      // std::cout<<std::endl;
      // std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddKuJ (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {

      // VALE, HO TINC BÉ, O COM A MÍNIM SEMBLANT EN EL MATLAB.
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      double dimension_double = double(dimension);

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      unsigned int voigtsize = 3;
      if (dimension == 3) 
         voigtsize = 6;

      Matrix Identity = ZeroMatrix( voigtsize, 1);

      for (unsigned int i = 0; i < dimension; i++) {
         Identity(i,0) = 1.0;
      }

      ConstitutiveMatrix = prod( ConstitutiveMatrix, (Identity) );
      ConstitutiveMatrix /= dimension_double;


      for ( unsigned int i = 0; i < voigtsize; i++)
      {
         ConstitutiveMatrix(i,0) += ( 2/dimension - 1.0 ) * rVariables.StressVector(i); 
      }

      double ElementJacobian = 0;
      for ( unsigned int i = 0; i <  number_of_nodes ; i++)
         ElementJacobian += GetGeometry()[i].GetSolutionStepValue( JACOBIAN ) * rVariables.N[i] ;


      ConstitutiveMatrix /= ElementJacobian;

      Matrix KuJ = prod ( trans( rVariables.B), ConstitutiveMatrix);
      Matrix SecondMatrix = ZeroMatrix( dimension*number_of_nodes, number_of_nodes);

      for (unsigned int i = 0; i < dimension*number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            SecondMatrix(i,j) = KuJ(i,0)*rVariables.N[j];
         }
      }

      SecondMatrix *= rIntegrationWeight;

      // ARA HE DE POSAR LA MATRIU AL SEU LLOC ( I AMB EL SEU SIGNE, NEGATIU??)
      MatrixType Kh=rLeftHandSideMatrix;
      unsigned int indexi = 0;
      unsigned int indexj = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int idim = 0; idim < dimension; idim++)
         {
            indexj = 0; 
            for (unsigned int j = 0; j < number_of_nodes; j++)
            {
               for (unsigned int jdim = 0; jdim < 1; jdim++)
               {
                  rLeftHandSideMatrix(indexi+i, indexj + 2*j+2) += SecondMatrix(indexi, indexj);
                  indexj++;
               }
            }
            indexi++;
         }
      }

      //std::cout << std::endl;
      //std::cout << " TRUE MATRIX " << rLeftHandSideMatrix - Kh << std::endl;
      //std::cout << std::endl;

      KRATOS_CATCH( "" )
   }



   //******************* Kuug ********************************************************
   //*********************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      int size = number_of_nodes * dimension;


      Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

      Matrix ReducedKg = prod( rVariables.DN_DX,  rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized

      Matrix Kuu = zero_matrix<double> (size);
      MathUtils<double>::ExpandAndAddReducedMatrix( Kuu, ReducedKg, dimension );

      // MatrixType Kh=rLeftHandSideMatrix;

      //assemble into rLeftHandSideMatrix the geometric uu contribution:
      unsigned int indexi = 0;
      unsigned int indexj = 0;
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

      // std::cout<<std::endl;
      // std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }



   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddKJu (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY



      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexp = dimension;


      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            int indexup = dimension*j + j;
            for (unsigned int k = 0; k < dimension; k++)
            {
               // Small Strains solid skeleton volume change
               rLeftHandSideMatrix(indexp, indexup + k) += rVariables.N[i] * rVariables.DN_DX( j, k) * rIntegrationWeight ;
               // LargeTerms solid skeleton volume change
            }
         }
         indexp += (dimension + 1);
      }


      //std::cout<<std::endl;
      //std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddKJJ (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix TotalF = prod( rVariables.F, rVariables.F0);

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;
      double consistent; 

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent = 1.0/12.0;
            if ( i == j)
               consistent *= 2.0;

            rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * rIntegrationWeight / rVariables.detFT;
            indexpj += (dimension+1);
         }

         indexpi += (dimension + 1);
      }

      KRATOS_CATCH( "" )
   }



   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::CalculateAndAddKJJStab (MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {

      KRATOS_TRY

      //repasar

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;
      double consistent = 1.0;

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

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=(-1.0)*AlphaStabilization;
            if(indexpi==indexpj)
               consistent=2.0*AlphaStabilization;

            rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);     //2D

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      // std::cout<<std::endl;
      // std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )

   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::GetHistoricalVariables( GeneralVariables& rVariables, const double& rPointNumber )
   {
      LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

      //Deformation Gradient F0
      rVariables.detF0 = mDeterminantF0[rPointNumber];
      rVariables.F0    = mDeformationGradientF0[rPointNumber];
   }

   //************************************CALCULATE VOLUME CHANGE*************************
   //************************************************************************************

   double& UpdatedLagrangianUJElement::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
   {
      KRATOS_TRY

      rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

      return rVolumeChange;

      KRATOS_CATCH( "" )
   }

   ////************************************************************************************
   ////************************************************************************************

   void UpdatedLagrangianUJElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      //create and initialize element variables:
      GeneralVariables Variables;
      this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

      //create constitutive law parameters:
      ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

      //set constitutive law flags:
      Flags &ConstitutiveLawOptions=Values.GetOptions();

      ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
      ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


      for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
      {

         //compute element kinematics B, F, DN_DX ...
         this->CalculateKinematics(Variables,PointNumber);

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


         //call the constitutive law to update material variables
         mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse(Values, Variables.StressMeasure);

         //call the constitutive law to finalize the solution step
         mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
               GetGeometry(),
               Variables.N,
               rCurrentProcessInfo );

         Variables.FT = ElementalFT;
         Variables.detFT = ElementalDetFT;


         //call the element internal variables update
         this->FinalizeStepVariables(Variables,PointNumber);
      }


      //PostProcess. Increment of Displacement
       Vector desp = ZeroVector(3);
         const unsigned int number_of_points = GetGeometry().PointsNumber();
         for (unsigned int i = 0; i < number_of_points; ++i)  {
         desp = ZeroVector(3);
         desp += GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
         desp -= GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT, 1);
         if ( GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT_INCR) )
         GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_INCR) = desp;
         } 

      mFinalizedStep = true;

      KRATOS_CATCH( "" )
   }
   ////************************************************************************************
   ////************************************************************************************

   ////************************************************************************************
   ////************************************************************************************

   ////************************************************************************************
   ////************************************************************************************

   void UpdatedLagrangianUJElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
         ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      //create and initialize element variables:
      GeneralVariables Variables;
      this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

      //create constitutive law parameters:
      ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

      //set constitutive law flags:
      Flags &ConstitutiveLawOptions=Values.GetOptions();

      ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
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

         //compute stresses and constitutive parameters
         mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, Variables.StressMeasure);

         Variables.FT = ElementalFT;
         Variables.detFT = ElementalDetFT;

         //some transformation of the configuration can be needed (UL element specially)
         this->TransformGeneralVariables(Variables,PointNumber);

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


         // std::cout<<" Element: "<<this->Id()<<std::endl;
         // unsigned int number_of_nodes = GetGeometry().PointsNumber();
         // for ( unsigned int i = 0; i < number_of_nodes; i++ )
         //   {
         //     array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
         //     array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
         //     array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
         //     array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
         //     std::cout<<" Previous  Position  node["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<std::endl;
         //   }
         // for ( unsigned int i = 0; i < number_of_nodes; i++ )
         //   {
         //     array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
         //     std::cout<<" Current  Position  node["<<GetGeometry()[i].Id()<<"]: "<<CurrentPosition<<std::endl;
         //   }
         // for ( unsigned int i = 0; i < number_of_nodes; i++ )
         //   {
         //     array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
         //     std::cout<<" Previous Displacement  node["<<GetGeometry()[i].Id()<<"]: "<<PreviousDisplacement<<std::endl;
         //   }

         // for ( unsigned int i = 0; i < number_of_nodes; i++ )
         //   {
         //     array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
         //     std::cout<<" Current  Displacement  node["<<GetGeometry()[i].Id()<<"]: "<<CurrentDisplacement<<std::endl;
         //   }
         // std::cout<<" Stress "<<Variables.StressVector<<std::endl;
         // std::cout<<" Strain "<<Variables.StrainVector<<std::endl;
         // std::cout<<" F  "<<Variables.F<<std::endl;
         // std::cout<<" F0 "<<Variables.F0<<std::endl;
         // std::cout<<" ConstitutiveMatrix "<<Variables.ConstitutiveMatrix<<std::endl;
         // std::cout<<" K "<<rLocalSystem.GetLeftHandSideMatrix()<<std::endl;
         // std::cout<<" f "<<rLocalSystem.GetRightHandSideVector()<<std::endl;



      }


      KRATOS_CATCH( "" )
   }


   void UpdatedLagrangianUJElement::ComputeConstitutiveVariables(  GeneralVariables& rVariables, Matrix& rFT, double& rDetFT)
      //void UpdatedLagrangianUJElement::ComputeConstitutiteVariables( const GeneralVariables& rVariables, Matrix& rFT, double& rDetFT)
   {
      KRATOS_TRY
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();

      rDetFT = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++) 
         rDetFT += GetGeometry()[i].GetSolutionStepValue( JACOBIAN ) * rVariables.N[i];

      rFT = rVariables.FT;

      rFT *=  pow( rDetFT/ rVariables.detFT, 1.0/double(dimension) );
      /*if ( this->Id() == 1) {
         std::cout << " IN THIS FUNCTION " << std::endl;
         std::cout << " FT " << rFT << std::endl;
         std::cout << " detFT " << rDetFT << std::endl;
         std::cout << " " << std::endl;
         std::cout << " VAARIAABLES " << std::endl;
         std::cout << " FT " << rVariables.FT << std::endl;
         std::cout << " DetFT " << rVariables.detFT << std::endl;
         std::cout << std::endl;
      }*/

      // COMPUTE THE EFFECT OF THE INTERPOLATION, LETS SEE
      std::vector< Matrix > EECCInverseDefGrad;
      ProcessInfo SomeProcessInfo;
      this->GetValueOnIntegrationPoints( INVERSE_DEFORMATION_GRADIENT, EECCInverseDefGrad, SomeProcessInfo);
      Matrix EECCInverseBig = EECCInverseDefGrad[0];
      Matrix EECCDefGradInverse = ZeroMatrix(dimension);
      for (unsigned int i = 0; i < dimension; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            EECCDefGradInverse(i,j) = EECCInverseBig(i,j);
         }
      }

      double det;
      Matrix EECCDefGrad;
      MathUtils<double>::InvertMatrix( EECCDefGradInverse, EECCDefGrad, det);

      double detF0 = 0;
      unsigned int step = 1;
      if ( mFinalizedStep ==  true) 
         step = 0;
         for ( unsigned int i = 0; i < number_of_nodes; i++)
         detF0 += GetGeometry()[i].GetSolutionStepValue( JACOBIAN, step ) * rVariables.N[i];

      Matrix F0 = rVariables.F0;
      F0 *= pow( detF0 / rVariables.detF0, 1.0/double(dimension) );

      Matrix F0Inverse;
      MathUtils<double>::InvertMatrix( F0, F0Inverse, det);

      Matrix Update = prod( F0Inverse, EECCDefGrad);

      if ( this->Id() == 0)
      {
         std::cout << " TRY TO SEE WHAT I DID " << std::endl;
         std::cout << " CONSTITUTIVE INVERSE " << EECCInverseDefGrad[0] << std::endl;
         std::cout << "  CONSTITUTIVE " << EECCDefGrad << std::endl;
         std::cout << std::endl;
         std::cout << " FINALIZED ?: " << mFinalizedStep << std::endl;
         std::cout << " NODAL " << detF0 << std::endl;
         std::cout << " PREVIOUS DISPL F " << rVariables.F0 << std::endl;
         std::cout << "   MAYBE " << rVariables.FT << std::endl;
         std::cout << " SO FINALLY F0 displ theta " << F0 << std::endl;
         std::cout << " and the Inverse is: " << F0Inverse << std::endl;
         std::cout << " UPDATE ?" << Update << std::endl;
         std::cout << " AND THE LAST ONE " << prod( rFT, Update) << std::endl;
         std::cout << std::endl;
      }

      // SO FINALLY I DO THAT
      rFT = prod( rFT, Update);


      KRATOS_CATCH( " " )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void UpdatedLagrangianUJElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }







}
