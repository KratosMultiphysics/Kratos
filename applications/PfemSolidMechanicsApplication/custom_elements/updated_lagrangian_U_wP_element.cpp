//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

//   Implementation of the Fluid Saturated porous media in a U-Pw Formulation
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_wP_element.hpp"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"
#include "custom_utilities/water_pressure_utilities.hpp" 

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   UpdatedLagrangianUwPElement::UpdatedLagrangianUwPElement()
      : UpdatedLagrangianElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPElement::UpdatedLagrangianUwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPElement::UpdatedLagrangianUwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUwPElement::UpdatedLagrangianUwPElement( UpdatedLagrangianUwPElement const& rOther)
      :UpdatedLagrangianElement(rOther)
       ,mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUwPElement&  UpdatedLagrangianUwPElement::operator=(UpdatedLagrangianUwPElement const& rOther)
   {
      UpdatedLagrangianElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPElement::~UpdatedLagrangianUwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void UpdatedLagrangianUwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ));
      }
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }
         else
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::GetValuesVector( Vector& rValues, int Step )
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
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
         }
         else
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
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

   void UpdatedLagrangianUwPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
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

   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUwPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-wP element type ", " UpdatedLagrangianUwPElement" )

            //verify that the variables are correctly initialized

            if ( PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }



   //**********************************GET DOUBLE VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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
      else{

         UpdatedLagrangianElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }

   }

   //**********************************GET VECTOR VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::GetValueOnIntegrationPoints( const Variable<Vector> & rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == DARCY_FLOW ) {

         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);

      }
      else{

         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }


   }

   //**********************************GET TENSOR VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo)
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
   void UpdatedLagrangianUwPElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
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
   void UpdatedLagrangianUwPElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable ==  DARCY_FLOW )
      {
         if ( rOutput.size() != mConstitutiveLawVector.size() ) {
            rOutput.resize( mConstitutiveLawVector.size() );
         }

         // CONSTITUTIVE PARAMETERS
         double Permeability = GetProperties()[PERMEABILITY];
         double WaterDensity = GetProperties()[DENSITY_WATER];

         // GEOMETRY PARAMETERS
         const unsigned int& integration_points_number = mConstitutiveLawVector.size();
         const unsigned int& dimension       = GetGeometry().WorkingSpaceDimension();
         const unsigned int& number_of_nodes = GetGeometry().size(); 

         // Get DN_DX
         ElementDataType Variables;
         this->InitializeElementData( Variables, rCurrentProcessInfo);

         Matrix K = ZeroMatrix( dimension, dimension);
         for (unsigned int i = 0; i < dimension; i++)
            K(i,i) = Permeability;  // this is only one of the two cases. 

         for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
         {
            this->CalculateKinematics(Variables, PointNumber);

            Vector GradP = ZeroVector( dimension );

            for (unsigned int i = 0; i < number_of_nodes; i++) {
               const double & rWaterPressure = GetGeometry()[i].FastGetSolutionStepValue( WATER_PRESSURE );
               for (unsigned int iDim = 0; iDim < dimension; iDim++) {
                  GradP(iDim) += Variables.DN_DX(i, iDim) * rWaterPressure; 
               }
            }

            // BTerm
            GradP(dimension-1) -= 10.0 * WaterDensity;

            // finally
            GradP  = prod( K, GradP);
            // Manual resize
            Vector ResizedVector = ZeroVector(3);
            for (unsigned int i = 0; i < dimension; i++)
               ResizedVector(i) = GradP(i);
            rOutput[PointNumber] = ResizedVector;
         }

      }
      else
      {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }
   }

   // *********************** Calculate Matrix On Integration Points *************************
   // ****************************************************************************************
   void UpdatedLagrangianUwPElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
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

   void UpdatedLagrangianUwPElement::InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
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

   void UpdatedLagrangianUwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
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


   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
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
      LocalSystemComponents BaseClassLocalSystem;
      Matrix BaseClassLeftHandSideMatrix = ZeroMatrix( dimension*number_of_nodes, dimension*number_of_nodes);
      BaseClassLocalSystem.SetLeftHandSideMatrix( BaseClassLeftHandSideMatrix );

      UpdatedLagrangianElement::CalculateAndAddLHS( BaseClassLocalSystem, rVariables, rIntegrationWeight);

      // Reshape the BaseClass LHS and Add the Hydro Part
      WaterPressureUtilities WaterUtility; 
      Matrix TotalF = prod ( rVariables.F, rVariables.F0);
      int number_of_variables = dimension + 1; 
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
      //HMVariables.CurrentRadius
      //HMVariables.ConstrainedModulus
      HMVariables.number_of_variables = number_of_variables;



      rLeftHandSideMatrix = WaterUtility.CalculateAndAddHydromechanicalLHS( HMVariables, rLeftHandSideMatrix, BaseClassLeftHandSideMatrix, rIntegrationWeight);


      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
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
      Vector BaseClassRightHandSideVector = ZeroVector ( dimension * number_of_nodes );
      BaseClassLocalSystem.SetRightHandSideVector(BaseClassRightHandSideVector );
      Vector VolumeForce = rVolumeForce; 
      VolumeForce *= 0.0;
      UpdatedLagrangianElement::CalculateAndAddRHS( BaseClassLocalSystem, rVariables, VolumeForce, rIntegrationWeight);


      // Reshape the BaseClass RHS and Add the Hydro Part
      WaterPressureUtilities WaterUtility;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);

      int number_of_variables = dimension+1; 


      // 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetDeformationGradient( TotalF);
      HMVariables.SetVolumeForce( rVolumeForce);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      //HMVariables.CurrentRadius
      //HMVariables.ConstrainedModulus
      HMVariables.number_of_variables = number_of_variables;


      rRightHandSideVector = WaterUtility.CalculateAndAddHydromechanicalRHS( HMVariables, rRightHandSideVector, BaseClassRightHandSideVector, rIntegrationWeight);


      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianElement )
   }

   void UpdatedLagrangianUwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianElement )
   }



}  // END KRATOS NAMESPACE
