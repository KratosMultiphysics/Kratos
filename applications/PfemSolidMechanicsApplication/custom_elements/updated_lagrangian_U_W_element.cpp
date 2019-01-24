//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                     PNavas $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                 October 2017 $
//   Revision:            $Revision:                     -0.1 $
//
//
//   Implementation of the Fluid Saturated porous media in a u-w Formulation

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_W_element.hpp"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"
#include "custom_utilities/water_pressure_utilities.hpp" 

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   UpdatedLagrangianUWElement::UpdatedLagrangianUWElement()
      : UpdatedLagrangianElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUWElement::UpdatedLagrangianUWElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUWElement::UpdatedLagrangianUWElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUWElement::UpdatedLagrangianUWElement( UpdatedLagrangianUWElement const& rOther)
      :UpdatedLagrangianElement(rOther)
       ,mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUWElement&  UpdatedLagrangianUWElement::operator=(UpdatedLagrangianUWElement const& rOther)
   {
      UpdatedLagrangianElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUWElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUWElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUWElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUWElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new UpdatedLagrangianUWElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUWElement::~UpdatedLagrangianUWElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void UpdatedLagrangianUWElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_DISPLACEMENT_Z ) );
      }
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = 2 * number_of_nodes * dimension;

      if ( rResult.size() != element_size )
         rResult.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int index = 2 * i * dimension;
         rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

         if( dimension == 3)
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_X ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_Y ).EquationId();
            rResult[index + 5] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_Z ).EquationId();
         } else {
            rResult[index + 2] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_X ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_Y ).EquationId();
         }


      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::GetValuesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = 2 * number_of_nodes * dimension;
      unsigned int dofs_per_node = dimension * 2;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dofs_per_node;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

         if ( dimension == 3 ) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_X, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Y, Step );
            rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Z, Step );
         } else {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_X, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Y, Step );
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = 2 * number_of_nodes * dimension;
      unsigned int dofs_per_node = dimension * 2;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dofs_per_node;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

         if ( dimension == 3 ) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_X, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Y, Step );
            rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Z, Step );
         } else {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_X, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Y, Step );
         }

      }
   }

   //*********************************ACCELERATION***************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = 2 * number_of_nodes * dimension;
      unsigned int dofs_per_node = dimension * 2;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dofs_per_node;

         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

         if ( dimension == 3 ) {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_X, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Y, Step );
            rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Z, Step );
         } else {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_X, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Y, Step );
         }

      }
   }

   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUWElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-wP element type ", " UpdatedLagrangianUWElement" )

            //verify that the variables are correctly initialized

            if ( PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  double WaterBulk = 1e+7;
      if ( GetProperties().Has(WATER_BULK_MODULUS)  ) {
         WaterBulk = GetProperties()[WATER_BULK_MODULUS];
      } else if ( rCurrentProcessInfo.Has(WATER_BULK_MODULUS) ) {
         WaterBulk = rCurrentProcessInfo[WATER_BULK_MODULUS];
      }
      GetProperties().SetValue(WATER_BULK_MODULUS, WaterBulk);

      double Permeability = 1e-5;
      if ( GetProperties().Has(PERMEABILITY)  ) {
         Permeability = GetProperties()[PERMEABILITY];
      } else if ( rCurrentProcessInfo.Has(PERMEABILITY) ) {
         Permeability = rCurrentProcessInfo[PERMEABILITY];
      }
      GetProperties().SetValue(PERMEABILITY, Permeability);

      double density = 0.0;
      if ( GetProperties().Has(DENSITY)  ) {
         density = GetProperties()[DENSITY];
      } else if ( rCurrentProcessInfo.Has(DENSITY) ) {
         density = rCurrentProcessInfo[DENSITY];
      }
      GetProperties().SetValue(DENSITY, density);

      double density_water = 0.0;
      if ( GetProperties().Has(DENSITY_WATER)  ) {
         density_water = GetProperties()[DENSITY_WATER];
      } else if ( rCurrentProcessInfo.Has(DENSITY_WATER) ) {
         density_water = rCurrentProcessInfo[DENSITY_WATER];
      }
      GetProperties().SetValue(DENSITY_WATER, density_water);

      double initial_porosity = 0.3;
      if ( GetProperties().Has(INITIAL_POROSITY) ) {
         initial_porosity = GetProperties()[INITIAL_POROSITY];
      } else if ( rCurrentProcessInfo.Has(INITIAL_POROSITY) ) {
         initial_porosity = rCurrentProcessInfo[INITIAL_POROSITY];
      }
      if ( initial_porosity < 1e-5)
         initial_porosity = 0.3;
      GetProperties().SetValue( INITIAL_POROSITY, initial_porosity);


      return correct;

      KRATOS_CATCH( "" );
   }



   //**********************************GET DOUBLE VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
         std::vector<double>& rValues,
         const ProcessInfo& rCurrentProcessInfo )
   {

      const unsigned int& integration_points_number = mConstitutiveLawVector.size();
      if ( rValues.size() != integration_points_number )
         rValues.resize( integration_points_number );

      if ( rVariable == POROSITY ||  rVariable == DENSITY )
      {
         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
      }
      else if ( rVariable == VOID_RATIO )
      {
         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
      }
      else if ( rVariable == WATER_PRESSURE )
      {
         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
      }
      else{

         UpdatedLagrangianElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }

   }

   //**********************************GET VECTOR VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::GetValueOnIntegrationPoints( const Variable<Vector> & rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
   {

      LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );



   }

   //**********************************GET TENSOR VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == TOTAL_CAUCHY_STRESS) 
      {
         CalculateOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);
      }
      else {

         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);

      }

   }


   // *********************** Calculate double On Integration Points *************************
   // ****************************************************************************************
   void UpdatedLagrangianUWElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {


      if ( rVariable == POROSITY ) {

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();
         const double InitialPorosity  = GetProperties()[INITIAL_POROSITY];

         if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize( mConstitutiveLawVector.size() );

         std::vector<double>  DetF0; 
         GetValueOnIntegrationPoints( DETERMINANT_F, DetF0, rCurrentProcessInfo);

         if (rOutput.size() != integration_points_number)
            rOutput.resize( integration_points_number) ;

         for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
         {
            rOutput[PointNumber] = 1.0 - (1.0 - InitialPorosity) / DetF0[PointNumber] ;
         }
      } else if (rVariable == WATER_PRESSURE ) {
         if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize( mConstitutiveLawVector.size() );

         ElementDataType Variables;
         this->InitializeElementData(Variables,rCurrentProcessInfo);
         const unsigned int PointNumber = 0;
         this->CalculateKinematics( Variables, PointNumber );

         double WaterPressure = this->CalculateGaussPointWaterPressure( Variables, WaterPressure);

         rOutput[PointNumber] = WaterPressure;

      }
      else if ( rVariable == VOID_RATIO) {

         GetValueOnIntegrationPoints( POROSITY, rOutput, rCurrentProcessInfo);

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();

         for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
         {
            rOutput[PointNumber] = rOutput[PointNumber] / (1.0 - rOutput[PointNumber]) ;
         }

      }
      else if ( rVariable == DENSITY )
      {
         const unsigned int& integration_points_number = mConstitutiveLawVector.size();

         if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize( mConstitutiveLawVector.size() );

         double DensityRigid, DensityWater, DensityMixtureInitial, DensityMixtureFinal, PorosityInitial, Porosity, Jacobian;

         DensityMixtureInitial = GetProperties()[DENSITY];
         DensityWater = GetProperties()[DENSITY_WATER];
         PorosityInitial = GetProperties()[INITIAL_POROSITY];

         DensityRigid = (DensityMixtureInitial-PorosityInitial*DensityWater) / (1.0 - PorosityInitial);

         Jacobian = mDeterminantF0[0];
         Porosity = 1.0 - ( 1.0 - PorosityInitial) / Jacobian; 

         DensityMixtureFinal = (1.0 - Porosity) * DensityRigid + Porosity * DensityWater;

         for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
         {
            rOutput[PointNumber] = DensityMixtureFinal;
         }

      }

      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
      }

   }

   // *********************** Calculate Vector On Integration Points *************************
   // ****************************************************************************************
   void UpdatedLagrangianUWElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
   }

   // *********************** Calculate Matrix On Integration Points *************************
   // ****************************************************************************************
   void UpdatedLagrangianUWElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == TOTAL_CAUCHY_STRESS)
      {
         if ( rOutput.size() != mConstitutiveLawVector.size() ) {
            rOutput.resize( mConstitutiveLawVector.size() );
         }
         // only for one gauss point element
         GetValueOnIntegrationPoints( CAUCHY_STRESS_TENSOR, rOutput, rCurrentProcessInfo);

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();
         const unsigned int number_of_nodes = GetGeometry().size();

         double GaussPointPressure = 0.0;
         unsigned int ThisDimension = rOutput[0].size1(); 
         Matrix Eye = ZeroMatrix(ThisDimension, ThisDimension);

         for (unsigned int i = 0; i < ThisDimension; ++i)
            Eye(i,i) = 1.0;

         for (unsigned int i = 0; i < number_of_nodes ; ++i)
            GaussPointPressure += GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE) ;

         GaussPointPressure /= double(number_of_nodes);

         for ( unsigned int PointNumber = 0; PointNumber <  integration_points_number; PointNumber++) {
            rOutput[PointNumber] = rOutput[PointNumber] + GaussPointPressure * Eye;
         }


      }
      else
      {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }
   }



   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      UpdatedLagrangianElement::InitializeElementData(rVariables,rCurrentProcessInfo);

      // SAVE THE TIME STEP, THAT WILL BE USED; BUT IS A BAD IDEA TO DO IT THIS WAY.
      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

      // KC permeability constitutive equation
      bool KozenyCarman = false; 
      if( GetProperties().Has(KOZENY_CARMAN) ){
         KozenyCarman = GetProperties()[KOZENY_CARMAN];
      }
      else if( rCurrentProcessInfo.Has(KOZENY_CARMAN) ){
         KozenyCarman = rCurrentProcessInfo[KOZENY_CARMAN];
      }
      GetProperties().SetValue(KOZENY_CARMAN, KozenyCarman);

      double initial_porosity  = 0.3; 
      if( GetProperties().Has(INITIAL_POROSITY) ){
         initial_porosity = GetProperties()[INITIAL_POROSITY];
      }
      else if( rCurrentProcessInfo.Has(INITIAL_POROSITY) ){
         KozenyCarman = rCurrentProcessInfo[INITIAL_POROSITY];
      }
      GetProperties().SetValue(INITIAL_POROSITY, initial_porosity);

   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize =  2 * number_of_nodes * dimension;

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

   void UpdatedLagrangianUWElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      KRATOS_TRY

      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      this->CalculateAndAddKww( rLeftHandSideMatrix, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //         Material stiffness matrix
   void UpdatedLagrangianUWElement::CalculateAndAddKuum( MatrixType & rLeftHandSide, ElementDataType& rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = dimension * 2;

      Matrix SmallMatrix = ZeroMatrix(number_of_nodes*dimension);

      noalias( SmallMatrix ) = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); 	

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            for (unsigned int k = 0; k < dimension; k++) {
               for (unsigned int l = 0; l < dimension; l++) {
                  rLeftHandSide( i*dofs_per_node + k, j*dofs_per_node + l) += SmallMatrix( i*dimension+k, j*dimension+l);
               }
            }
         }
      }



      KRATOS_CATCH("")
   }

   //************************************************************************************
   //      Water material Matrix
   void UpdatedLagrangianUWElement::CalculateAndAddKww( MatrixType & rLeftHandSide, ElementDataType& rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY


      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = dimension * 2;

      unsigned int shear_size=1;
      if (dimension == 3) shear_size=3;
      unsigned int mat_B2_size = 2*dimension+shear_size;

      Matrix B2 = ZeroMatrix(mat_B2_size, number_of_nodes*(2*dimension));
      CalculateB2Matrix( B2, rVariables.DN_DX );   //No valdria para 3D PNA

      Matrix Q = ZeroMatrix(mat_B2_size, mat_B2_size);
      double Bulk = GetProperties()[WATER_BULK_MODULUS]; 
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);
      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0;
      Bulk /= porosity;


      for (unsigned int i = 0; i < dimension; i++)
         for (unsigned int j = 0; j < dimension; j++){
            Q(i,j) = Bulk;
            Q(dimension+shear_size+i,j) = Bulk;
            Q(i,dimension+shear_size+j) = Bulk;
            Q(dimension+shear_size+i,dimension+shear_size+j) = Bulk;
         }

      Matrix SmallMatrix( number_of_nodes*dofs_per_node, number_of_nodes*dofs_per_node);
      noalias( SmallMatrix ) = prod( trans(B2), rIntegrationWeight * Matrix( prod( Q, B2) ) );
      rLeftHandSide += SmallMatrix;

      //std::cout << "Kel_con agua = " << rLeftHandSide << std::endl;

      KRATOS_CATCH("")
   }


   //*********************************************************************************
   //   Calculate some sort of B matrix

   void UpdatedLagrangianUWElement::CalculateB2Matrix( Matrix & rB2, const Matrix & rDN_DX)
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         if ( dimension == 2){
            rB2( 0, i*2*dimension) = rDN_DX(i,0);
            rB2( 1, i*2*dimension+1) = rDN_DX(i,1);
            rB2( 2, i*2*dimension+1) = rDN_DX(i,0);
            rB2( 2, i*2*dimension) = rDN_DX(i,1);
            rB2( dimension+1, i*2*dimension+2) = rDN_DX(i,0);
            rB2( dimension+2, i*2*dimension+3) = rDN_DX(i,1);
         }
      }



      KRATOS_CATCH("")

   }

   // ***********************************************************************************
   // ***********************************************************************************
   void UpdatedLagrangianUWElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      KRATOS_TRY

      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;


      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 



      this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight); 

      this->CalculateAndAddExternalWaterForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight); 

      this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      this->CalculateAndAddWaterPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);



      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }


   // **************************************************************************
   // Calculate and Add volumetric loads
   void UpdatedLagrangianUWElement::CalculateAndAddExternalForces( VectorType & rRightHandSideVector, ElementDataType & rVariables, 
         Vector & rVolumeForce, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = dimension * 2;

      rVolumeForce *= rVariables.detF0; 
      double density_mixture0 = GetProperties().GetValue(DENSITY);
      if ( density_mixture0 > 0) {
         rVolumeForce /= density_mixture0; 
      }
      else {
         return; 
      }

      double density_water =GetProperties().GetValue(DENSITY_WATER);
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);

      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 
      double density_solid = (density_mixture0 - porosity0*density_water) / ( 1.0 - porosity0);
      double density_mixture = ( 1.0 - porosity) * density_solid + porosity * density_water;


      const VectorType & rN = rVariables.N;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            rRightHandSideVector(i*dofs_per_node + iDim) += rIntegrationWeight * density_mixture *  rN(i) * rVolumeForce(iDim);
         }
      }

      rVolumeForce /= rVariables.detF0; 
      rVolumeForce *=density_mixture0;
      return;


      KRATOS_CATCH("")
   }

   // **************************************************************************
   // Calculate and Add volumetric loads
   void UpdatedLagrangianUWElement::CalculateAndAddExternalWaterForces( VectorType & rRightHandSideVector, ElementDataType & rVariables, 
         Vector & rVolumeForce, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = dimension * 2;

      rVolumeForce *= rVariables.detF0; 
      double density_mixture0 = GetProperties().GetValue(DENSITY);
      if ( density_mixture0 > 0) {
         rVolumeForce /= density_mixture0; 
      }
      else {
         return; 
      }

      double density_water =GetProperties().GetValue(DENSITY_WATER);

      const VectorType & rN = rVariables.N;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            rRightHandSideVector(i*dofs_per_node+dimension + iDim) += rIntegrationWeight * density_water *  rN(i) * rVolumeForce(iDim);
         }
      }

      rVolumeForce /= rVariables.detF0;
      rVolumeForce *=density_mixture0;
      return;


      KRATOS_CATCH("")
   }

   // **********************************************************************************
   //         CalculateAndAddInternalForces
   void UpdatedLagrangianUWElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = dimension * 2;

      double WaterPressure = this->CalculateGaussPointWaterPressure( rVariables, WaterPressure);


      Vector StressVector = rVariables.StressVector;
      for (unsigned int i = 0; i < dimension; i++)
         StressVector(i) -= WaterPressure;

      VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );


      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            rRightHandSideVector(i*dofs_per_node + j) -= InternalForces(i*dimension + j);
         }
      }

      KRATOS_CATCH( "" )
   }


   // **********************************************************************************
   //         CalculateAndAddWater-like-InternalForces
   void UpdatedLagrangianUWElement::CalculateAndAddWaterPressureForces( VectorType & rRightHandSideVector, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = dimension * 2;

      double WaterPressure = this->CalculateGaussPointWaterPressure( rVariables, WaterPressure);


      unsigned int voigt_size = 3;
      if ( dimension == 3)
         voigt_size = 6;

      Vector WaterStressVector = ZeroVector( voigt_size);
      for (unsigned int i = 0; i < dimension; i++)
         WaterStressVector(i) -= WaterPressure;

      //std::cout << "B = " << rVariables.B << std::endl;
      VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), WaterStressVector );


      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            rRightHandSideVector(i*dofs_per_node + j+ dimension) -= InternalForces(i*dimension + j);
         }
      }

      KRATOS_CATCH( "" )
   }

   // *********************************************************************************
   //        Calculate the Water Pressure at integration point
   double & UpdatedLagrangianUWElement::CalculateGaussPointWaterPressure( ElementDataType & rVariables, double & rWaterPressure)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double divU = 0;
      double divW = 0;

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            const array_1d<double, 3> & rCurrentDisplacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3> & rCurrentWaterDisplacement = GetGeometry()[i].GetSolutionStepValue(WATER_DISPLACEMENT);
            divU += rCurrentDisplacement[j]      * rVariables.DN_DX(i,j);
            divW += rCurrentWaterDisplacement[j] * rVariables.DN_DX(i,j);
         }
      }

      double Bulk = GetProperties()[WATER_BULK_MODULUS]; 

      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);
      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0;
      Bulk /= porosity;

      rWaterPressure = - Bulk * ( divU + divW);

      return rWaterPressure;

      KRATOS_CATCH("")
   }


   // *********************************************************************************
   //         Calculate the Mass matrix
   void UpdatedLagrangianUWElement::CalculateMassMatrix( MatrixType & rMassMatrix, ProcessInfo & rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension;
      unsigned int MatSize = 2 * number_of_nodes * dimension;

      if ( rMassMatrix.size1() != MatSize )
         rMassMatrix.resize( MatSize, MatSize, false );

      rMassMatrix = ZeroMatrix( MatSize, MatSize );


      //reading integration points
      IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::GI_GAUSS_2; //GeometryData::GI_GAUSS_1;

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      double density_mixture0 = GetProperties()[DENSITY];
      double WaterDensity =GetProperties().GetValue(DENSITY_WATER);
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);

      double porosity = 1.0 - (1.0-porosity0) / Variables.detF0; 
      double density_solid = (density_mixture0 - porosity0*WaterDensity) / ( 1.0 - porosity0);
      double CurrentDensity = ( 1.0 - porosity) * density_solid + porosity * WaterDensity;


      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
         //compute element kinematics
         this->CalculateKinematics( Variables, PointNumber );

         //getting informations for integration
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;

         IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );  // multiplies by thickness



         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               for ( unsigned int k = 0; k < dimension; k++ )
               {
                  rMassMatrix( dofs_per_node*i+k, dofs_per_node*j +k  )                     += Variables.N[i] * Variables.N[j] * CurrentDensity * IntegrationWeight;
                  rMassMatrix( dofs_per_node*i+dimension+k, dofs_per_node*j +k  )           += Variables.N[i] * Variables.N[j] * WaterDensity   * IntegrationWeight;
                  rMassMatrix( dofs_per_node*i+k, dofs_per_node*j+dimension +k  )           += Variables.N[i] * Variables.N[j] * WaterDensity   * IntegrationWeight;
                  rMassMatrix( dofs_per_node*i+dimension+k, dofs_per_node*j+dimension +k  ) += Variables.N[i] * Variables.N[j] * WaterDensity   * IntegrationWeight / porosity;
               }
            }
         }

      }


      KRATOS_CATCH("")
   }


   // *********************************************************************************
   //         Calculate the Damping matrix
   void UpdatedLagrangianUWElement::CalculateDampingMatrix( MatrixType & rDampingMatrix, ProcessInfo & rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension;
      unsigned int MatSize = 2 * number_of_nodes * dimension;

      if ( rDampingMatrix.size1() != MatSize )
         rDampingMatrix.resize( MatSize, MatSize, false );

      rDampingMatrix = ZeroMatrix( MatSize, MatSize );


      //reading integration points
      IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::GI_GAUSS_2; //GeometryData::GI_GAUSS_1;

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);



      double CurrentPermeability = GetProperties()[PERMEABILITY]; 

      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
         //compute element kinematics
         this->CalculateKinematics( Variables, PointNumber );

         //std::cout << "N = " << Variables.N << std::endl;

         //getting informations for integration
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;

         IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );  // multiplies by thickness


         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               for ( unsigned int k = 0; k < dimension; k++ )
               {
                  rDampingMatrix( dofs_per_node*i+dimension+k, dofs_per_node*j+dimension +k  ) += Variables.N[i] * Variables.N[j] * IntegrationWeight / CurrentPermeability;
               }
            }
         }

      }

      KRATOS_CATCH("")
   }
   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUWElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianElement )
   }

   void UpdatedLagrangianUWElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianElement )
   }



}  // END KRATOS NAMESPACE
