//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                     PNavas $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                 October 2017 $
//   Revision:            $Revision:                     -0.1 $
//
//
//   Implementation of the Fluid Saturated porous media in a u-w-pw Formulation

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_U_J_W_wP_element.hpp"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   AxisymUpdatedLagrangianUJWwPElement::AxisymUpdatedLagrangianUJWwPElement()
      : AxisymUpdatedLagrangianUJElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPElement::AxisymUpdatedLagrangianUJWwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : AxisymUpdatedLagrangianUJElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPElement::AxisymUpdatedLagrangianUJWwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : AxisymUpdatedLagrangianUJElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPElement::AxisymUpdatedLagrangianUJWwPElement( AxisymUpdatedLagrangianUJWwPElement const& rOther)
      :AxisymUpdatedLagrangianUJElement(rOther)
       ,mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPElement&  AxisymUpdatedLagrangianUJWwPElement::operator=(AxisymUpdatedLagrangianUJWwPElement const& rOther)
   {
      AxisymUpdatedLagrangianUJElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJWwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUJWwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJWwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUJWwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new AxisymUpdatedLagrangianUJWwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPElement::~AxisymUpdatedLagrangianUJWwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void AxisymUpdatedLagrangianUJWwPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
   {
      rElementalDofList.resize( 0 );

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( JACOBIAN) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_DISPLACEMENT_Y ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
	}
      
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJWwPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      unsigned int element_size          = number_of_nodes * 6;
      unsigned int dofs_per_node = 6;

      if ( rResult.size() != element_size )
         rResult.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int index =  i  * dofs_per_node;
         rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_X ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( WATER_DISPLACEMENT_Y ).EquationId();
            rResult[index + 5] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJWwPElement::GetValuesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      unsigned int element_size          = number_of_nodes * 6;
      unsigned int dofs_per_node = 6;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dofs_per_node;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_X, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_DISPLACEMENT_Y, Step );
            rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJWwPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      unsigned int element_size          = number_of_nodes * 6;
      unsigned int dofs_per_node = 6;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dofs_per_node;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
            rValues[index + 2] = 0.0;
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_X, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_VELOCITY_Y, Step );
            rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_VELOCITY, Step );
      }
   }

   //*********************************ACCELERATION***************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJWwPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      unsigned int element_size          = number_of_nodes * 6;
      unsigned int dofs_per_node = 6;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dofs_per_node;

         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[index + 2] = 0.0;
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_X, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_ACCELERATION_Y, Step );
            rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_ACCELERATION, Step );
      }
   }

   //************************************************************************************
   //************************************************************************************

   int  AxisymUpdatedLagrangianUJWwPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = AxisymUpdatedLagrangianUJElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-wP element type ", " AxisymUpdatedLagrangianUJWwPElement" )

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

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJWwPElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      AxisymUpdatedLagrangianUJElement::InitializeElementData(rVariables,rCurrentProcessInfo);

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

   void AxisymUpdatedLagrangianUJWwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize =  number_of_nodes * ( 2*dimension + 2);

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

   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int dofs_per_node = 2*dimension + 2;

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

      // Assemble the Kuj Left Hand Side
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            for (unsigned int iDim = 0; iDim < dimension+1; iDim++) {
               for (unsigned int jDim = 0; jDim < dimension+1; jDim++) {
                  rLeftHandSideMatrix( i*dofs_per_node + iDim, j*dofs_per_node+jDim) += LocalLeftHandSideMatrix(i*(dimension+1) + iDim, j*(dimension+1) + jDim );
               }
            }
         }
      }

      this->CalculateAndAddKWwP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      this->CalculateAndAddKUwP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      this->CalculateAndAddKPPStab( rLeftHandSideMatrix, rVariables, IntegrationWeight);
      
      this->CalculateAndAddHighOrderKPP( rLeftHandSideMatrix, rVariables, IntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }

   // **********************************************************************************
   //          Matrix that may appear due to the stabilization matrix 
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddKPPStab( MatrixType & rLeftHandSideMatrix,ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // **********************************************************************************
   //          Matrix that may appear due to the high order terms
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddHighOrderKPP( MatrixType & rLeftHandSideMatrix,ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //         Matrix due to the the water pressure contribution to the internal forces   
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddKUwP( MatrixType & rLeftHandSide, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

      Matrix Q = ZeroMatrix(number_of_nodes*dimension, number_of_nodes);

      unsigned int voigtSize = 4;
      Matrix m = ZeroMatrix( voigtSize, 1);

      for ( unsigned int i = 0; i < 3; i++)
         m(i,0) = 1.0;
      Matrix partial =  prod( trans( rVariables.B), m );

      for (unsigned int i = 0; i < dimension*number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++ ) {
            Q(i,j) = partial(i,0) * rVariables.N[j] * rIntegrationWeight;
         }
      }

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               rLeftHandSide( i*dofs_per_node + iDim, (j+1)*dofs_per_node-1) -= Q( i*dimension+iDim, j);
            }
         }
      }


      KRATOS_CATCH("")
   }


   //************************************************************************************
   //      Water material Matrix
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddKWwP( MatrixType & rLeftHandSide, ElementDataType& rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

      Matrix SmallMatrix = ZeroMatrix(number_of_nodes*dimension, number_of_nodes);

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               SmallMatrix( i*dimension+iDim, j ) = rVariables.N[i] * rVariables.DN_DX(j,iDim) * rIntegrationWeight;
               //if (iDim==0)
               //   SmallMatrix( i*dimension+iDim, j ) += rVariables.N[i] * rVariables.N(j)/rVariables.CurrentRadius * rIntegrationWeight;
            }
         }
      }


      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            for (unsigned int iDim = 0; iDim < dimension; iDim++) {
               rLeftHandSide( i*dofs_per_node +(dimension+1)+iDim, (j+1)*dofs_per_node-1) += SmallMatrix( i*dimension+iDim, j);
            }
         }
      }



      KRATOS_CATCH("")
   }



   // ***********************************************************************************
   // ***********************************************************************************
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      KRATOS_TRY

      // define some variables
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int dofs_per_node = 2 * (dimension+1);

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

      // Assemble the "UJ" RHS
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension + 1; iDim++) {
            rRightHandSideVector( i*dofs_per_node + iDim) += BaseClassRightHandSideVector( i*(dimension+1) + iDim);
         }
      }


      this->CalculateAndAddExternalForcesUJWwP( rRightHandSideVector, rVariables, rVolumeForce, IntegrationWeight); 

      this->CalculateAndAddInternalWaterForces( rRightHandSideVector, rVariables, IntegrationWeight);


      // Calculate the linear momentum of the fluid
      CalculateAndAddFluidLinearMomentum( rRightHandSideVector, rVariables, IntegrationWeight);

      // Calculate mass balance equation
      CalculateAndAddMassBalanceEquation( rRightHandSideVector, rVariables, IntegrationWeight);

      this->CalculateAndAddStabilizationRHS( rRightHandSideVector, rVariables, IntegrationWeight);
      this->CalculateAndAddHighOrderRHS( rRightHandSideVector, rVariables, IntegrationWeight);



      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }


   // *********************************************************************************
   //    part of the RHS that goes directly to the RHS
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddStabilizationRHS( VectorType & rRightHandSideVector, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // *********************************************************************************
   //    High order part of the RHS that goes directly to the RHS
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddHighOrderRHS( VectorType & rRightHandSideVector, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // **********************************************************************************
   //    mass balance equation of the mixture (aka: Darcy's Law ) 
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddMassBalanceEquation( VectorType & rRightHandSideVector, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY
      // a convective term may go here. not coded yet. 
      KRATOS_CATCH("")
   }

   // **********************************************************************************
   //    linear momentum balance equation of the fluid phase (aka: Darcy's Law ) 
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddFluidLinearMomentum( VectorType & rRightHandSideVector, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = 2*dimension + 2;


      Vector GradP = ZeroVector(dimension);
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         const double & rWaterPressure = GetGeometry()[i].FastGetSolutionStepValue(WATER_PRESSURE);
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            GradP(iDim ) += rVariables.DN_DX(i,iDim) * rWaterPressure;
         }
      }


      Vector SmallRHS = ZeroVector(dimension*number_of_nodes);
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            SmallRHS(i*dimension+iDim) = rVariables.N(i)*GradP(iDim) * rIntegrationWeight;
         }
      }

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            rRightHandSideVector(i*dofs_per_node + (dimension+1)+iDim ) -= SmallRHS(i*dimension+iDim);
         }
      }

      // LMV. a term due to the water pressure weight is missing. Please coded ASAP!!

      KRATOS_CATCH("")
   }

   // **************************************************************************
   // Calculate and Add volumetric loads
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddExternalForcesUJWwP( VectorType & rRightHandSideVector, ElementDataType & rVariables, 
         Vector & rVolumeForce, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node = 2*dimension + 2;

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

   // **********************************************************************************
   //         CalculateAndAddInternalForces
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddInternalWaterForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int dofs_per_node =  2 * dimension + 2;


      double WaterPressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         WaterPressure += GetGeometry()[i].FastGetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];

      unsigned int voigt_size = 4;
      Vector StressVector = ZeroVector(voigt_size);
      for (unsigned int i = 0; i < 3; i++)
         StressVector(i) -= WaterPressure;

      VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );


      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            rRightHandSideVector(i*dofs_per_node + j) -= InternalForces(i*dimension + j);
         }
      }

      KRATOS_CATCH( "" )
   }



   // *********************************************************************************
   //         Calculate the Mass matrix
   void AxisymUpdatedLagrangianUJWwPElement::CalculateMassMatrix( MatrixType & rMassMatrix, ProcessInfo & rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension + 2;
      unsigned int MatSize =  number_of_nodes * ( 2 * dimension + 2);

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
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ *2.0 * 3.141592654 * Variables.CurrentRadius;

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               for ( unsigned int k = 0; k < dimension; k++ )
               {
                  rMassMatrix( dofs_per_node*i+k, dofs_per_node*j +k  )                     += Variables.N[i] * Variables.N[j] * CurrentDensity * IntegrationWeight;
                  rMassMatrix( dofs_per_node*i+(dimension+1)+k, dofs_per_node*j +k  )           += Variables.N[i] * Variables.N[j] * WaterDensity   * IntegrationWeight;
                  rMassMatrix( dofs_per_node*i+k, dofs_per_node*j+(dimension+1) +k  )           += Variables.N[i] * Variables.N[j] * WaterDensity   * IntegrationWeight;
                  rMassMatrix( dofs_per_node*i+(dimension+1)+k, dofs_per_node*j+(dimension+1) +k  ) += Variables.N[i] * Variables.N[j] * WaterDensity   * IntegrationWeight / porosity;
               }
            }
         }

         this->CalculateAndAddMassStabilizationMatrix( rMassMatrix, Variables, IntegrationWeight);
      }


      KRATOS_CATCH("")
   }

   // ********************************************************************************
   //      part of the mass matrix that steams from the stabilization factor
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddMassStabilizationMatrix( MatrixType & rMassMatrix, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // ********************************************************************************
   //      part of the damping matrix coming from the high order terms
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddHighOrderDampingMatrix( MatrixType & rDampingMatrix, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   // *********************************************************************************
   //         Calculate the Damping matrix
   void AxisymUpdatedLagrangianUJWwPElement::CalculateDampingMatrix( MatrixType & rDampingMatrix, ProcessInfo & rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension + 2;
      unsigned int MatSize = number_of_nodes * ( 2*dimension + 2 );

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
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ *2.0 * 3.141592654 * Variables.CurrentRadius;


         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               for ( unsigned int k = 0; k < dimension; k++ )
               {
                  rDampingMatrix( dofs_per_node*i+(dimension+1)+k, dofs_per_node*j+(dimension+1) +k  ) += Variables.N[i] * Variables.N[j] * IntegrationWeight / CurrentPermeability;
               }
            }
         }


         // Q Matrix //

         Matrix Q = ZeroMatrix( number_of_nodes, dimension*number_of_nodes);
         unsigned int voigtSize = 4;
         Matrix m = ZeroMatrix( 1, voigtSize);
         for ( unsigned int i = 0; i < 3; i++)
            m(0,i) = 1.0;
         Matrix partial =  prod( m, Variables.B );
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < dimension*number_of_nodes; j++ ) {
               Q(i,j) = Variables.N[i] * partial(0,j) * IntegrationWeight;
            }
         }

         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               for (unsigned int jDim = 0; jDim < dimension; jDim++) {
                  rDampingMatrix( (i+1)*dofs_per_node -1, j*dofs_per_node + jDim ) += Q(i, j*dimension+jDim);
                  rDampingMatrix( (i+1)*dofs_per_node -1, j*dofs_per_node + jDim + (dimension+1) ) += Q(i, j*dimension+jDim);
               }
            }
         }

         Matrix Mass = ZeroMatrix( number_of_nodes, number_of_nodes);
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               Mass(i,j) += Variables.N[i]*Variables.N[j] * IntegrationWeight;
            }
         }
         const double  rWaterBulk = GetProperties()[WATER_BULK_MODULUS];
         Mass /= rWaterBulk;

         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               rDampingMatrix( (i+1)*dofs_per_node-1, (j+1)*dofs_per_node-1) += Mass(i,j);
            }
         }

         this->CalculateAndAddDampingStabilizationMatrix(rDampingMatrix, Variables, IntegrationWeight);
         this->CalculateAndAddHighOrderDampingMatrix(rDampingMatrix, Variables, IntegrationWeight);


      } // end point
      KRATOS_CATCH("")
   }


   // ********************************************************************************
   //      part of the damping matrix that steams from the stabilization factor
   void AxisymUpdatedLagrangianUJWwPElement::CalculateAndAddDampingStabilizationMatrix( MatrixType & rDampingMatrix, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension + 2;
      
      
      const double & rStabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];
      if ( ( fabs(rStabilizationFactor) > 1.0e-6) && dimension==2)  {

         double StabFactor = CalculateStabilizationFactor( rVariables, StabFactor);


         Matrix SmallMatrix = ZeroMatrix(number_of_nodes, number_of_nodes);

         double consistent;
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               consistent = -1.0 * StabFactor / 18.0;
               if ( i == j)
                  consistent = 2.0 * StabFactor / 18.0;
               SmallMatrix(i,j) += consistent * rIntegrationWeight ;
            }
         }

         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               rDampingMatrix( (i+1)*dofs_per_node-1, (j+1)*dofs_per_node-1) += SmallMatrix(i,j);
            }
         }
      }

      KRATOS_CATCH("")
   }
   //************************************************************************************
   //************************************************************************************

   double & AxisymUpdatedLagrangianUJWwPElement::CalculateStabilizationFactor( ElementDataType & rVariables, double & rStabFactor)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ElementSize = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         double aux = 0;
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            aux += rVariables.DN_DX(i, iDim);
         }
         ElementSize += fabs( aux); 
      }
      ElementSize *= sqrt( double(dimension) );
      ElementSize = 4.0/ ElementSize; 

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


      double StabilizationFactor = GetProperties().GetValue( STABILIZATION_FACTOR_WP);


      //rStabFactor = 2.0 / ConstrainedModulus - 12.0 * rPermeability * mTimeStep / pow(ElementSize, 2); 
      //rStabFactor = 2.0 / ConstrainedModulus*(1-rPermeability*DryDensity/mTimeStep/porosity) - 12.0 * rPermeability * mTimeStep / pow(ElementSize, 2)*(1-porosity)/porosity; 



      rStabFactor = 2.0 / ConstrainedModulus; // - 12.0 * rPermeability * mTimeStep / pow(ElementSize, 2); 
      rStabFactor *=  StabilizationFactor;

      if ( rStabFactor < 0)
	      rStabFactor = 0;

      return rStabFactor; 

      KRATOS_CATCH("")
   }

   void AxisymUpdatedLagrangianUJWwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUJElement )
   }

   void AxisymUpdatedLagrangianUJWwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUJElement )
   }



}  // END KRATOS NAMESPACE
