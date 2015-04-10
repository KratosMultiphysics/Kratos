//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/large_displacement_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"
#include "custom_elements/spatial_lagrangian_U_wP_element.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   SpatialLagrangianUwPElement::SpatialLagrangianUwPElement()
      : LargeDisplacementElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPElement::SpatialLagrangianUwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : LargeDisplacementElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPElement::SpatialLagrangianUwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : LargeDisplacementElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SpatialLagrangianUwPElement::SpatialLagrangianUwPElement( SpatialLagrangianUwPElement const& rOther)
      :LargeDisplacementElement(rOther)
       ,mDeformationGradientF0(rOther.mDeformationGradientF0)
       ,mDeterminantF0(rOther.mDeterminantF0)
       ,mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   SpatialLagrangianUwPElement&  SpatialLagrangianUwPElement::operator=(SpatialLagrangianUwPElement const& rOther)
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

   Element::Pointer SpatialLagrangianUwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new SpatialLagrangianUwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      SpatialLagrangianUwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
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

      return Element::Pointer( new SpatialLagrangianUwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPElement::~SpatialLagrangianUwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void SpatialLagrangianUwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
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

   void SpatialLagrangianUwPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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

   void SpatialLagrangianUwPElement::GetValuesVector( Vector& rValues, int Step )
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

   void SpatialLagrangianUwPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
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

   void SpatialLagrangianUwPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
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

   int  SpatialLagrangianUwPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

         int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_ERROR( std::logic_error, "constitutive law is not compatible with the U-wP element type ", " SpatialLagrangianUwPElement" )

            //verify that the variables are correctly initialized

            if ( PRESSURE.Key() == 0 )
               KRATOS_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }

   //*********************************SET DOUBLE VALUE***********************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
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


   void SpatialLagrangianUwPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
         std::vector<double>& rValues,
         const ProcessInfo& rCurrentProcessInfo )
   {

      if (rVariable == DETERMINANT_F){

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();

         if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

         for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
         {
            rValues[PointNumber] = mDeterminantF0[PointNumber];
         }

      }
      else if ( rVariable == POROSITY)
      {
         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
      }
      else{

         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }

   }



   void SpatialLagrangianUwPElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
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
      else if ( rVariable == POROSITY ) {

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();
         const double InitialPorosity  = GetProperties()[INITIAL_POROSITY];

         std::vector<double>  DetF0; 
         CalculateOnIntegrationPoints( DETERMINANT_F, DetF0, rCurrentProcessInfo);

         if (rOutput.size() != integration_points_number)
            rOutput.resize( integration_points_number) ;

         for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
         {
            rOutput[PointNumber] = 1.0 - (1.0 - InitialPorosity) / DetF0[PointNumber] ;
         }
         //std::cout << " ELEMENT: " <<  this->Id() << " POROSITY " << rOutput[0] << std::endl;
         //std::cout << " INITIAL PO " << InitialPorosity << " and " << 1.0 - (1.0-InitialPorosity)/DetF0[0] << std::endl;
      }
      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }

   }



   void SpatialLagrangianUwPElement::GetValueOnIntegrationPoints( const Variable<Vector> & rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == DARCY_FLOW ) {

         double ScalingConstant;
         double Permeability; double WaterBulk; double DeltaTime;

         GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();


         const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
         Vector b = ZeroVector(dimension);
         double WaterDensity = GetProperties()[DENSITY_WATER];
         b(dimension-1) = -10.0*WaterDensity;

         if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

         const GeometryType& rGeom = this->GetGeometry();
         //const unsigned int NumberOfNodes = rGeom.PointsNumber(); 

         GeometryType::ShapeFunctionsGradientsType DN_DX;
         Vector GaussWeights;
         Vector DetJ;
         rGeom.ShapeFunctionsIntegrationPointsGradients( DN_DX, DetJ, GeometryData::GI_GAUSS_1);
         //    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
         //    {
         unsigned int PointNumber = 0;
         const Matrix& ThisPointDN_DX = DN_DX[PointNumber];
         //        std::cout << " I AM HERE " << DN_DX << std::endl;
         //        std::cout << std::endl;
         Vector Pressures(3);
         for (unsigned int i = 0; i < 3 ; ++i)
            Pressures(i) = GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE);
         Vector Darcy = prod ( trans(ThisPointDN_DX), Pressures);
         Vector D1 = ZeroVector(3);

         for (unsigned int i  = 0; i < 2; ++i)
            D1(i) = (Darcy(i) + b(i));
         rValues[PointNumber] = D1;

         for (unsigned int i = 0; i < integration_points_number; i++) {
            rValues[i] = D1 * Permeability;
         }
         //       std::cout << " PRESSURES " << Pressures << std::endl;
         //       std::cout << std::endl;
         //       std::cout << prod( ThisPointDN_DX, Pressures) << std::endl;
         //       std::cout << prod( trans(Pressures), trans(ThisPointDN_DX)) <<std::endl;
         //       std::cout << prod( trans(ThisPointDN_DX), Pressures) << std::endl;
         //       std::cout << std::endl;
         //    }
      }
      else{

         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

      }


   }
   //**********************************GET TENSOR VALUE**********************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo)
   {
      if (rVariable == TOTAL_CAUCHY_STRESS) {

         GetValueOnIntegrationPoints( CAUCHY_STRESS_TENSOR, rValue, rCurrentProcessInfo);
         //std::cout << " I AM SUPPOSED TO BE HERE " << rValue << std::endl;

         const unsigned int& integration_points_number = mConstitutiveLawVector.size();
         //const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

         double GaussPointPressure = 0.0;
         unsigned int ThisDimension = rValue[0].size1(); 
         Matrix Eye = ZeroMatrix(ThisDimension);

         for (unsigned int i = 0; i < ThisDimension; ++i)
            Eye(i,i) = 1.0;

         for (unsigned int i = 0; i < 3 ; ++i)
            GaussPointPressure += GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE);

         GaussPointPressure /= 3.0;

         for ( unsigned int PointNumber = 0; PointNumber <  integration_points_number; PointNumber++) {
            rValue[PointNumber] = rValue[PointNumber] + GaussPointPressure * Eye;
         }

      }
      else {

         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);

      }


   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************
   void SpatialLagrangianUwPElement::Initialize()
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


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      LargeDisplacementElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

      //Calculate Delta Position
      rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

      //set variables including all integration points values

      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

      // SAVE THE TIME STEP, THAT WILL BE USED; BUT IS A BAD IDEA TO DO IT THIS WAY.
      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

      //mCompressibleWater = false;

   }

   ////************************************************************************************
   ////************************************************************************************

   void SpatialLagrangianUwPElement::FinalizeStepVariables( GeneralVariables & rVariables, const double& rPointNumber )
   { 
      //update internal (historical) variables
      mDeterminantF0[rPointNumber]         = rVariables.detF0 ;
      mDeformationGradientF0[rPointNumber] = rVariables.F0;
   }


   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
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


   void SpatialLagrangianUwPElement::CalculateDeformationMatrix(Matrix& rB,
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

         KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

      }

      KRATOS_CATCH( "" )
   }



   //*************************COMPUTE DEFORMATION GRADIENT*******************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::CalculateDeformationGradient(const Matrix& rDN_DX,
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

         KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

      }

      KRATOS_CATCH( "" )
   }

   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************
   //************************************************************************************
   //************************************************************************************

   //************************************************************************************
   void SpatialLagrangianUwPElement::SetGeneralVariables(GeneralVariables& rVariables,
         ConstitutiveLaw::Parameters& rValues,
         const int & rPointNumber)
   {
      LargeDisplacementElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

      //Set extra options for the contitutive law
      Flags &ConstitutiveLawOptions=rValues.GetOptions();
      ConstitutiveLawOptions.Set(ConstitutiveLaw::FINAL_CONFIGURATION);

   }

   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************
   void SpatialLagrangianUwPElement::CalculateKinematics(GeneralVariables& rVariables,
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

      //Step domain size
      rVariables.DomainSize = rVariables.detJ;

      //Compute cartesian derivatives [dN/dx_n]
      noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

      //Current Deformation Gradient [dx_n+1/dx_n]
      //this->CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition);

      //Deformation Gradient F [dx_n+1/dx_n] to be updated
      noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

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

   void SpatialLagrangianUwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
   {


      /* std::cout << " I HAVE IT ALL WRONG: detF0 " << 1.0 - rVariables.detF0 << " detF " << 1.0 - rVariables.detF << std::endl;
         std::cout << " PART 2 F  " << rVariables.F << std::endl;
         std::cout << " PART 3 F0 " << rVariables.F0 << std::endl;
         std::cout << std::endl;
       */


      // EM FALTA UN CACHO
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
      CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpu to the rLefsHandSideMatrix
      CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpp to the rLefsHandSideMatrix
      CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpp Stab to the rLefsHandSideMatrix
      CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;




      //KRATOS_WATCH( rLeftHandSideMatrix )
   }

   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
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
      CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }


   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
         GeneralVariables& rVariables,
         Vector& rVolumeForce,
         double& rIntegrationWeight)

   {
      KRATOS_TRY
         unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // VectorType Fh=rRightHandSideVector;

      double DomainSize = (rVariables.DomainSize / rVariables.detJ );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int indexup = dimension * i + i;
         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] += rIntegrationWeight * rVariables.N[i] * rVolumeForce[j] * DomainSize;
         }

      }

      // std::cout<<std::endl;
      // std::cout<<" Fext "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )
   }


   //************************** INTERNAL FORCES    *******************************
   //************************************** Idem but with Total Stress ***********

   void SpatialLagrangianUwPElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // VectorType Fh=rRightHandSideVector;

      Vector TotalStressVector = rVariables.StressVector;

      double Pressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; ++i) {
         Pressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];
      }

      for (unsigned int i = 0; i < dimension; ++i) 
         TotalStressVector(i) += Pressure;


      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), TotalStressVector );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexup = dimension * i + i;
         unsigned int indexu  = dimension * i;

         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
         }
      }

      // std::cout<<std::endl;
      // std::cout<<" Stress "<<rVariables.StressVector<<std::endl;
      // std::cout<<" Fint "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )
   }

   //******************************** PRESSURE FORCES  **********************************
   //***************************************** aka: Mass conservation equation **********

   void SpatialLagrangianUwPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      Vector b = ZeroVector(dimension);
      double WaterDensity = GetProperties()[DENSITY_WATER];
      b(dimension-1) = -10.0*WaterDensity;
      // WHY??

      Matrix K;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);
      this->GetPermeabilityTensor( Permeability, TotalF, K);

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

         indexp += (dimension + 1);

      }

    // std::cout<<std::endl;
    // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
    // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;
      KRATOS_CATCH( "" )

   }



   //****************** STABILIZATION *********************************************************
   //************************* defined in the Stab element ************************************

   void SpatialLagrangianUwPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
   }


   //******** Kuu Material************************************************************
   //***************** It includes the pw geometric stiffness ************************

   void SpatialLagrangianUwPElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         //assemble into rk the material uu contribution:
         const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // FIRST: Ctotal = Cmaterial + pJ(1*1 - 2I4)
      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;

      double Pressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; ++i)
         Pressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];

      unsigned int voigtsize = 3;
      if (dimension == 3) 
         voigtsize = 6;

      Matrix FourthOrderIdentity = ZeroMatrix(voigtsize);
      Matrix MatrixProduct = ZeroMatrix(voigtsize);
      for (unsigned int i = 0; i < dimension ; ++i) {
         FourthOrderIdentity(i,i) = 1.0;
         for (unsigned int j = 0; j < dimension; ++j) {
            MatrixProduct(i,j) = 1.0;
         }
      }
      for (unsigned int i = dimension; i < voigtsize; ++i)
         FourthOrderIdentity(i,i) = 0.5;


      //    ConstitutiveMatrix += Pressure * rVariables.detF0 * (MatrixProduct - 2.0*FourthOrderIdentity); // CREC QUE EL DET NO FA FALTA, ja és cauchy
      ConstitutiveMatrix += Pressure *  (MatrixProduct - 2.0*FourthOrderIdentity); 

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




   //******************* Kuug ********************************************************
   //*********************************************************************************

   void SpatialLagrangianUwPElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      int size = number_of_nodes * dimension;


      double Pressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; ++i)
         Pressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );
      Pressure /= double(number_of_nodes);


      Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

      for (unsigned int i = 0; i < dimension; ++i)
         StressTensor(i,i) += Pressure; //TotalStress

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

   void SpatialLagrangianUwPElement::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //MatrixType Kh=rLeftHandSideMatrix;
      //contributions to stiffness matrix calculated on the reference configuration
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexp  = dimension;
         unsigned int indexup = dimension * i + i;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;
            }
            indexp += (dimension + 1);
         }
      }

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
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
      this->GetPermeabilityTensor( Permeability, TotalF, K);

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexp = dimension;

      // try to compute Grad( u^n+1 - u^n)
      Matrix Velocity_DX = ZeroMatrix(dimension);
      Vector Pressure_DX = ZeroVector(dimension);
      Matrix Identity = ZeroMatrix(dimension);

      /// COMPUTE SOME GRADIENTS TO THE WaterPressure and Velocity
      // since there is only one gauss point there is no need to put the grad(u) and thinks like that as u_i grad(N),... etc.
      for (unsigned int k = 0; k < number_of_nodes; k++)
      {
         const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[k].FastGetSolutionStepValue( DISPLACEMENT );
         const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[k].FastGetSolutionStepValue( DISPLACEMENT , 1 );
         array_1d<double, 3 > DeltaDisplacement            = CurrentDisplacement-PreviousDisplacement;
         const double & WaterPressure = GetGeometry()[k].FastGetSolutionStepValue( WATER_PRESSURE );
         for (unsigned int j = 0; j < dimension; j++)
         {
            for (unsigned int i = 0; i < dimension; i++)
            {
               Velocity_DX(i,j) += DeltaDisplacement[i] * rVariables.DN_DX(k,j);

            }
            Pressure_DX(j) += WaterPressure * rVariables.DN_DX(k,j);
         }
      }

      for (unsigned int k = 0; k < dimension; k++)
      { 
         Identity(k,k) = 1.0;
      }


      Matrix FTotal = prod( rVariables.F, rVariables.F0);
      double LDTerm;

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            int indexup = dimension*j + j;
            for (unsigned int k = 0; k < dimension; k++)
            {
               // Small Strains solid skeleton volume change
               rLeftHandSideMatrix(indexp, indexup + k) += rVariables.N[i] * rVariables.DN_DX( j, k) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               // LargeTerms solid skeleton volume change
               for (unsigned int l = 0; l < dimension; l++)
               {
                  rLeftHandSideMatrix(indexp, indexup + k) -= rVariables.N[i] * Velocity_DX(l,k) * rVariables.DN_DX(j,l) * rIntegrationWeight * ScalingConstant / rVariables.detF0;

                  for (unsigned n = 0; n < dimension; n++)
                  {
                     for (unsigned t = 0; t < dimension ; t++)
                     {
                        // Permeability LD Terms (non-linear terms)
                        LDTerm =this->GetPermeabilityLDTerm( K, TotalF, l,n,k,t);
                        LDTerm = - LDTerm +  Identity(k,l)*K(t,n) + Identity(k,n)*K(t,l);
                        rLeftHandSideMatrix(indexp, indexup + k ) += rVariables.DN_DX(i,l) * Pressure_DX(n) * ( LDTerm ) * rVariables.DN_DX(j,t) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
                     }
                  }
               }
            }
         }
         indexp += (dimension + 1);
      }


      // std::cout<<std::endl;
      // std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      Matrix K;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);
      this->GetPermeabilityTensor( Permeability, TotalF, K);

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(WaterBulk)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * ScalingConstant / rVariables.detF0;
            for ( unsigned int p = 0; p < dimension; ++p ) {
               for ( unsigned int q = 0; q < dimension; ++q ) {
                  rLeftHandSideMatrix(indexpi, indexpj) -= DeltaTime * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q) * rIntegrationWeight* ScalingConstant / rVariables.detF0;
               }
            }

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      KRATOS_CATCH( "" )
   }



   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
   }


   //************************************************************************************
   // **** TO BE DESTROYED BECAUSE IT DOES NOT MAKE ANY SENSE ******
   //************************************************************************************
   void SpatialLagrangianUwPElement::GetConstants(double& rScalingConstant, double& rWaterBulk, double& rDeltaTime, double& rPermeability)
   {
      //double ScalingConstant = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

      rScalingConstant = 0.10;
      //rScalingConstant = 10.0;
      //rWaterBulk = 1.0e+8*ScalingConstant;
      rDeltaTime  = mTimeStep;

      rPermeability = GetProperties()[PERMEABILITY];
      rWaterBulk = GetProperties()[WATER_BULK_MODULUS];

   }


   // function to get the Eulerian permeability assuming a constant eulerian permeability
   void SpatialLagrangianUwPElement::GetPermeabilityTensor( const double& rPermeability, const Matrix& rF, Matrix& rPermeabilityTensor)
   {

      unsigned int thisSize = rF.size1();
      rPermeabilityTensor = ZeroMatrix( thisSize);
      for (unsigned int i = 0; i < thisSize; ++i)
      {
         rPermeabilityTensor(i,i) = rPermeability;
      }
      return;

      rPermeabilityTensor = prod(rPermeabilityTensor, trans( rF) );
      rPermeabilityTensor = prod( rF, rPermeabilityTensor);
      rPermeabilityTensor /= MathUtils<double>::Det(rF);

   }

   double SpatialLagrangianUwPElement::GetPermeabilityLDTerm( const Matrix& rPermeability, const Matrix& rF, const int i, const int j, const int k, const int l)
   {
      return 0.0;

      unsigned int thisSize = rF.size1();
      Matrix Identity = ZeroMatrix(thisSize);
      for (unsigned int h = 0; h < thisSize; h++)
         Identity(h,h) = 1.0;

      double LDTerm;
      LDTerm = Identity(i,k)*rPermeability(j,l) + Identity(i,l)*rPermeability(j,k);
      LDTerm = -rPermeability(i,j)*Identity(k,l);
      return LDTerm;

   }
   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void SpatialLagrangianUwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }







}
