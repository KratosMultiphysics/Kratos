//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_P_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"

//
// U-J-P FORMULATION OF CHAPTER 5, VOL 2 OF ZIEN, ED 6.
// THE RESIDUAL WORKS
//
// APPROXIMATED VERSION WITHOUT THE TERM NodalJ/ElemJ
//
// (even in plane strain it requeres the constitutive tensor in 3D (6X6) and the stress (1X6) 
//

//
namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUJPElement::UpdatedLagrangianUJPElement()
      : UpdatedLagrangianUJElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJPElement::UpdatedLagrangianUJPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUJElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJPElement::UpdatedLagrangianUJPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUJElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUJPElement::UpdatedLagrangianUJPElement( UpdatedLagrangianUJPElement const& rOther)
      :UpdatedLagrangianUJElement(rOther)
       //,mDeformationGradientF0(rOther.mDeformationGradientF0)
       //,mDeterminantF0(rOther.mDeterminantF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUJPElement&  UpdatedLagrangianUJPElement::operator=(UpdatedLagrangianUJPElement const& rOther)
   {
      UpdatedLagrangianUJElement::operator=(rOther);

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

   Element::Pointer UpdatedLagrangianUJPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUJPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUJPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUJPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJPElement::~UpdatedLagrangianUJPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   void UpdatedLagrangianUJPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
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
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
      }
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = number_of_nodes * dimension + 2*number_of_nodes; //LL

      if ( rResult.size() != element_size )
         rResult.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int index = i * dimension + 2*i;
         rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

         if( dimension == 3)
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
         }
         else
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
         }

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::GetValuesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension + 2*i;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
         }
         else
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension + 2*i;
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

   void UpdatedLagrangianUJPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * dimension + 2*i;
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

   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUJPElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = UpdatedLagrangianUJElement::Check(rCurrentProcessInfo);

      //verify that the variables are correctly initialized

      if ( JACOBIAN.Key() == 0 )
         KRATOS_THROW_ERROR( std::invalid_argument, "JACOBIAN has Key zero! (check if the application is correctly registered", "" )
            if ( PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }


   //**********************************GET DOUBLE VALUE**********************************
   //************************************************************************************


   void UpdatedLagrangianUJPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
         std::vector<double>& rValues,
         const ProcessInfo& rCurrentProcessInfo )
   {

      UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

   }

   void UpdatedLagrangianUJPElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {

      KRATOS_TRY


      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rOutput.size() != integration_points_number )
         rOutput.resize( integration_points_number );

      //if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) 
      // SINCE IN THIS VARIABLE I PUT THE "CONSTITUTIVE STRESS" It DOESN?T MATTER

      UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);

      KRATOS_CATCH("")
   }


   void UpdatedLagrangianUJPElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      
      if ( rVariable == CAUCHY_STRESS_TENSOR ) {
         const unsigned int number_of_nodes = GetGeometry().PointsNumber();

         //create and initialize element variables:
         ElementDataType Variables;
         this->InitializeElementData(Variables, rCurrentProcessInfo);

         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions = Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

         // for integration points
         for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {

            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //to take in account previous step writing
            if( this->Is(SolidElement::FINALIZED_STEP) ){
               this->GetHistoricalVariables(Variables,PointNumber);
            }		

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

            double NodalPressure = 0;
            for (unsigned int i = 0; i < number_of_nodes; i++) {
               NodalPressure += GetGeometry()[i].GetSolutionStepValue(PRESSURE) * Variables.N[i];
            }

         // OBS, now changing Variables I change Values because they are pointers ( I hope);
            double ElementalDetFT = Variables.detH;
            Matrix ElementalFT = Variables.H;

            // AND NOW IN THE OTHER WAY
            Matrix m; double d;
            ComputeConstitutiveVariables( Variables, m, d);

            Variables.H = m;
            Variables.detH = d;
            Values.SetDeformationGradientF( Variables.H );
            Values.SetDeterminantF( Variables.detH );


            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);

            // T1
            Variables.H  = ElementalFT;
            Variables.detH = ElementalDetFT;

            Vector StressVector = Variables.StressVector;
            double ElementalPressure = 0;
            for (unsigned int i = 0; i < 3; i++)
               ElementalPressure += StressVector(i);
            ElementalPressure /= 3.0;

            //NodalPressure *= ( Variables.detF/ NodalJacobian);
            for (unsigned int i = 0; i < 3; i++)
               StressVector(i) += ( NodalPressure - ElementalPressure);

            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector);
         }

      }
      else {
         UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }
   }

   void UpdatedLagrangianUJPElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {

      UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);

   }



   void UpdatedLagrangianUJPElement::GetValueOnIntegrationPoints( const Variable<Vector> & rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
   {

      UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
   }

   //**********************************GET TENSOR VALUE**********************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo)
   {
      if ( rVariable == CAUCHY_STRESS_TENSOR)
      {
         CalculateOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);
      }
      else {

         UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);
      }

   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************
   void UpdatedLagrangianUJPElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      UpdatedLagrangianUJElement::InitializeElementData(rVariables,rCurrentProcessInfo);

      rVariables.StressVector.resize(6);

      rVariables.ConstitutiveMatrix.resize(6,6);


   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
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

   void UpdatedLagrangianUJPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {


      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // operation performed: add Km to the rLefsHandSideMatrix
      UJPElementData  ElementVariables; 
      CalculateThisElementData( ElementVariables, rVariables);

      //respect to the current configuration n+1
      CalculateAndAddKuumElemUJP( rLeftHandSideMatrix, rVariables, ElementVariables,  rIntegrationWeight );

      // operation performed: add Kg to the rLefsHandSideMatrix
      CalculateAndAddKuugElemUJP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      CalculateAndAddKuJElemUJP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kup to the rLefsHandSideMatrix
      CalculateAndAddKup( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );


      CalculateAndAddKJuElemUJP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      CalculateAndAddKJJElemUJP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      CalculateAndAddKJp( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpu to the rLefsHandSideMatrix
      CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      CalculateAndAddKpJ( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpp to the rLefsHandSideMatrix
      CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpp Stab to the rLefsHandSideMatrix
      CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );
      CalculateAndAddKJJStabElemUJP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      //std::cout << " SYSTEMMATRIX " << rLeftHandSideMatrix << std::endl;

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;


   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      if (this->Id() == 0 ) {
         std::cout << " FT " << rVariables.detH << std::endl;
         std::cout << " FF " << rVariables.detF << std::endl;
         std::cout << " F0 " << rVariables.detF0 << std::endl;
         std::cout << " " << std::endl;
      }

      // EM FALTA UN CACHO
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      UJPElementData  ElementVariables; 
      CalculateThisElementData( ElementVariables, rVariables);


      // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
      CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

      // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
      CalculateAndAddInternalForcesElemUJP( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      CalculateAndAddJacobianForcesElemUJP( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
      CalculateAndAddPressureForces( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      CalculateAndAddStabilizedJacobianElemUJP( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;



      //KRATOS_WATCH( rRightHandSideVector )
   }



   //************************** INTERNAL FORCES    *******************************
   //************************************** Idem but with Total Stress ***********

   void UpdatedLagrangianUJPElement::CalculateAndAddInternalForcesElemUJP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         UJPElementData& rElementVariables, 
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      VectorType Fh=rRightHandSideVector;

      Vector StressVector = rElementVariables.StressVector;

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexup = dimension * i + 2*i;
         unsigned int indexu  = dimension * i;

         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
         }
      }

      KRATOS_CATCH( "" ) 
   }

   //******************************** PRESSURE FORCES  **********************************
   //************************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddJacobianForcesElemUJP( VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         UJPElementData& rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension ;

      VectorType Fh=rRightHandSideVector;

      double JacobianElement = rVariables.detH;

      double consistent = 1.0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent = 1.0/12.0;
            if ( i==j)
               consistent *= 2.0;
            const double& JacobianNodal = GetGeometry()[j].GetSolutionStepValue(JACOBIAN) ;

            rRightHandSideVector[indexp] -=  consistent * JacobianNodal * rIntegrationWeight / rVariables.detH ;

         }

         rRightHandSideVector[indexp] += rVariables.N[i] * JacobianElement * rIntegrationWeight / rVariables.detH;


         indexp += (dimension + 2);
      }

      KRATOS_CATCH( "" )

   }

   //******************************** PRESSURE FORCES  **********************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         UJPElementData& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension + 1;

      VectorType Fh=rRightHandSideVector;

      double ElementalMeanStress = rElementVariables.ElementalMeanStress; 
      double consistent = 1.0;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent = 1.0/12.0;
            if ( i == j)
               consistent *= 2.0;
            const double& Pressure = GetGeometry()[j].GetSolutionStepValue(PRESSURE) ;
            rRightHandSideVector[indexp] -=   consistent * Pressure * rIntegrationWeight / rVariables.detH ;

         }

         rRightHandSideVector[indexp] += rVariables.N[i] * ElementalMeanStress * rIntegrationWeight / rVariables.detH;

         indexp += (dimension + 2);

      }

      KRATOS_CATCH( "" )
   }



   //****************** STABILIZATION *********************************************************
   //************************* defined in the Stab element ************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddStabilizedJacobianElemUJP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      //VectorType Fh = rRightHandSideVector;

      double AlphaStabilization  = 4.0; 
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
      AlphaStabilization *= StabilizationFactor; 


      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

      AlphaStabilization *= BulkModulus;  // TIMES THE BULK MODULUS BECAUSE I HAVE ALL THE EQUATION MULTIPLIED BY THE BULK MODULUS
      if ( YoungModulus < 0.0001)
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0;
         AlphaStabilization *= mElementStabilizationNumber; 

      }

      double consistent = 1.0;

      unsigned int indexp = dimension ;
      for ( unsigned int i = 0; i < number_of_nodes; i++)
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++)
         {
            consistent = (-1.0)*AlphaStabilization;
            if ( i==j )
               consistent = 2.0*AlphaStabilization;

            const double& DetFNodal = GetGeometry()[j].GetSolutionStepValue(JACOBIAN);
            rRightHandSideVector[indexp] -= consistent * DetFNodal * rIntegrationWeight / (rVariables.detF0 / rVariables.detF);
         }
         indexp += (dimension + 2);
      }

      KRATOS_CATCH( "" )
   }


   // **************************** ADD STABILIZED PRESSURE *************************
   // ********************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      //VectorType Fh = rRightHandSideVector;

      double AlphaStabilization = 4.0;
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_P];
      AlphaStabilization *= StabilizationFactor;

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2.0*(1.0+PoissonCoefficient));
      double Bulk = YoungModulus/ ( 3.0 * ( 1.0 - 2.0*PoissonCoefficient) );

      AlphaStabilization = (4.0*StabilizationFactor/(18.0*LameMu)) * Bulk;  // TIMES THE BULK MODULUS
      if ( YoungModulus < 0.0001)
         return;

      double consistent = 1.0;

      unsigned int indexp = dimension + 1;
      for ( unsigned int i = 0; i < number_of_nodes; i++)
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++)
         {
            consistent = (-1.0)*AlphaStabilization;
            if ( i==j )
               consistent = 2.0*AlphaStabilization;

            const double& Pressure = GetGeometry()[j].GetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] -= consistent * Pressure * rIntegrationWeight / rVariables.detF0;
         }
         indexp += (dimension + 2);
      }

      KRATOS_CATCH( "" )

   }


   //******** Kuu Material************************************************************
   //***************** It includes the pw geometric stiffness ************************

   void UpdatedLagrangianUJPElement::CalculateAndAddKuumElemUJP(MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData &  rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      Matrix ConstMatrixSmall= ZeroMatrix(rElementVariables.voigtsize,rElementVariables.voigtsize) ;

      // Definition of Smallt Cons Matrix
      if ( rElementVariables.voigtsize == 6) {
         ConstMatrixSmall = rVariables.ConstitutiveMatrix; 
      }
      else {
         for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j< 3; j++) {
               unsigned int indexi = i;
               unsigned int indexj = j;
               if ( i== 2) indexi += 1;
               if ( j== 2) indexj += 1;
               ConstMatrixSmall(i, j)  = rVariables.ConstitutiveMatrix(indexi,indexj) ;
            }
         }
      }

      // 2. Definition of Deviatoric Alpha and multiply
      Matrix DeviatoricAlpha = ZeroMatrix ( rElementVariables.voigtsize);
      for (unsigned int i = 0; i < rElementVariables.voigtsize; i++)
         DeviatoricAlpha(i,i) = 1.0;
      for (unsigned int i = 0; i < dimension; i++) {
         for (unsigned int j = 0; j< dimension; j++) {
            DeviatoricAlpha(i,j) -= rElementVariables.Alpha;
         }
      }

      ConstMatrixSmall = prod( ConstMatrixSmall, DeviatoricAlpha);


      // 3. Some Geometric-like terms ( obs, copied from the other side);
      Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector);
      Matrix Identity = ZeroMatrix(3,3);
      for (unsigned int i = 0; i < dimension; i++)
         Identity(i,i) = 1.0;

      Matrix MoreTerms = ZeroMatrix(6,6);
      for ( unsigned int k = 0; k < dimension; k++) {
         for (unsigned int l = 0; l < dimension; l++) {
            for (unsigned int p = 0; p < dimension ; p++) {
               for (unsigned int q = 0; q < dimension; q++) {
                  // THIS IS NOT THE WAY TO PERFORM A 4Or mutliplication and compression, but,....
                  unsigned int auxk, auxl;
                  if ( k < l) { auxk = k; auxl = l;  } 
                  else  { auxk = l; auxl = k; }
                  unsigned int indexi;
                  if ( auxk == auxl) {
                     indexi = auxk;
                  }  else {
                     if ( auxk == 1) {
                        if ( auxl == 2) {
                           indexi = 3;
                        } else {
                           indexi = 4;
                        }
                     } else {
                        indexi = 5;
                     }
                  }
                  unsigned int auxp, auxq;
                  if ( p < q) { auxp = p;  auxq = q;}
                  else  {  auxp = q; auxq = p; }
                  unsigned int indexj;
                  if ( auxp == auxq) {
                     indexj = auxp;
                  } else {
                     if ( auxp == 1) {
                        if ( auxq == 2) {
                           indexj = 3;
                        } else {
                           indexj = 4;
                        }
                     } else {
                        indexj = 5;
                     }
                  }
                  double voigtnumber = 1.0;
                  if ( indexi > 2)
                     voigtnumber *= 0.5;
                  if (indexj > 2)
                     voigtnumber *= 0.5;

                  MoreTerms(indexi, indexj) += voigtnumber * ( Identity(k,p) * StressTensor(l,q) +  Identity(l,p) * StressTensor(k,q) - 2.0 * rElementVariables.Alpha * StressTensor(k,l)*Identity(p,q) );
               }
            }
         }
      }



      if ( this->Id() == 0) {
         std::cout << " CONST MATRIX SMALL " << ConstMatrixSmall << std::endl;
         std::cout << " MORE TERMS " << MoreTerms << std::endl;
      }

      // 4. Put it in the BIG MATRIX ( add out of plane terms) and add the ExtraTerms ( "geometric-like") 
      Matrix ConstMatrixBig = ZeroMatrix(6,6);
      if ( rElementVariables.voigtsize == 6) {
         ConstMatrixBig = ConstMatrixSmall + MoreTerms; 
      } 
      else {
         for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j< 3; j++) {
               unsigned int indexi = i;
               unsigned int indexj = j;
               if ( i== 2) indexi += 1;
               if ( j== 2) indexj += 1;
               ConstMatrixBig(indexi, indexj)  = ConstMatrixSmall(i,j) + MoreTerms(indexi, indexj) ;
            }
         }
         if ( this->Id() == 0) {
            std::cout << " I ADD THE OUT OF PLANE "<< ConstMatrixBig  << std::endl; 
         }

         // AND NOW LETS ADD THE OUT OF PLANE
         Matrix Aux1 = ZeroMatrix(1,3);
         for (unsigned int i = 0; i < 2; i++) {
            unsigned int indexi = i;
            if ( i == 2) indexi += 1;
            Aux1(0,i) = rVariables.ConstitutiveMatrix(2, indexi);
         }
         if (this->Id() == 0 ) {
            std::cout << " PREIVIOUS TO MULT " << Aux1 << std::endl;
         }
         Aux1 = prod ( Aux1, DeviatoricAlpha);

         if (this->Id() == 0) {
            std::cout << " AUX 1"  << Aux1 << std::endl; 
         }
         for (unsigned int i = 0; i < 2; i++) {
            unsigned int indexi = i;
            ConstMatrixBig(2, indexi) = Aux1(0, i);
         }
      }

      if ( this->Id() == 0) {
         std::cout << " AFTER THE OUT OF PLANE (i.e. breve derivative) " << ConstMatrixBig << std::endl;
      }

      //5. Multiply it by the Deviatoric Beta Matrix 
      Matrix DeviatoricBeta = ZeroMatrix(6,6);
      for (unsigned int i = 0; i < 6; i++) 
         DeviatoricBeta(i,i) = 1.0;
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3 ; j++) {
            DeviatoricBeta(i,j) -= rElementVariables.Beta; 
         }
      }
      if ( this->Id() == 0 ) {
         std::cout << " PREVIOS MULTIPLICATION " << std::endl;
         std::cout << " CONST MATRIX " << ConstMatrixBig << std::endl;
         std::cout << " DEVIATORIC BETA " << DeviatoricBeta << std::endl;
      }

      ConstMatrixBig = prod( DeviatoricBeta, ConstMatrixBig);

      // 6. Put the matrix in the small matrix
      if ( rElementVariables.voigtsize == 6) {
         ConstMatrixSmall = ConstMatrixBig;
      }
      else {
         for (unsigned int i = 0; i < rElementVariables.voigtsize; i++) {
            for (unsigned int j = 0; j < rElementVariables.voigtsize; j++) {
               unsigned int indexi = i;
               unsigned int indexj = j;
               if ( i == 2) indexi += 1;
               if ( j == 2) indexj += 1;
               ConstMatrixSmall(i,j) = ConstMatrixBig(indexi, indexj); 
            }
         }
      }



      // 7. Put the matrix in the MATRIX
      //contributions to stiffness matrix calculated on the reference config
      Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( ConstMatrixSmall, rVariables.B ) ) );  

      if (this->Id() == 0) {
         std::cout << " DISPLACEMENT DEVIATORIC DERIVATIVE " << ConstMatrixSmall << std::endl;
         std::cout << " THIS DIMENSIONS " << Kuu.size1() << " and " << Kuu.size2() << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
      }

      // MatrixType Kh=rLeftHandSideMatrix;

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
                  rLeftHandSideMatrix(indexi+2*i,indexj+2*j)+=Kuu(indexi,indexj);
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

   void UpdatedLagrangianUJPElement::CalculateAndAddKuugElemUJP(MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)

   {

      // Matrix from the paper, it is not the appropiate way to do it, but.....
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();



      Matrix ThisMatrix = ZeroMatrix(6,6);
      Matrix ThisMatrixSize = ZeroMatrix(rElementVariables.voigtsize,rElementVariables.voigtsize);

      ///// MALAMENT. ES LA TENSIÓ GUAY; NO UNA QUALSEVOL. SAPS??
      // MIRA AL PAPER.
      Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector);

      for (unsigned int i = 0; i < 3; i++)
         StressTensor(i,i) += (  rElementVariables.NodalMeanStress - rElementVariables.ElementalMeanStress);

      Matrix Identity = ZeroMatrix(3,3);
      for (unsigned int i = 0; i < dimension; i++)
         Identity(i,i) = 1.0;

      for ( unsigned int k = 0; k < dimension; k++) {
         for (unsigned int l = 0; l < dimension; l++) {
            for (unsigned int p = 0; p < dimension ; p++) {
               for (unsigned int q = 0; q < dimension; q++) {
                  // THIS IS NOT THE WAY TO PERFORM A 4Or mutliplication and compression, but,....
                  unsigned int auxk, auxl;
                  if ( k < l) { 
                     auxk = k; auxl = l;
                  } else  {
                     auxk = l; auxl = k;
                  }
                  unsigned int indexi;
                  if ( auxk == auxl) {
                     indexi = auxk;
                  }
                  else {
                     if ( auxk == 1) {
                        if ( auxl == 2) {
                           indexi = 3;
                        } else {
                           indexi = 4;
                        }
                     } else {
                        indexi = 5;
                     }
                  }
                  unsigned int auxp, auxq;
                  if ( p < q) { 
                     auxp = p;  auxq = q;
                  } else  {
                     auxp = q; auxq = p;
                  }
                  unsigned int indexj;
                  if ( auxp == auxq) {
                     indexj = auxp;
                  }
                  else {
                     if ( auxp == 1) {
                        if ( auxq == 2) {
                           indexj = 3;
                        } else {
                           indexj = 4;
                        }
                     } else {
                        indexj = 5;
                     }
                  }
                  double voigtnumber = 1.0;
                  if ( indexi > 2)
                     voigtnumber *= 0.5;
                  if (indexj > 2)
                     voigtnumber *= 0.5;

                  ThisMatrix(indexi, indexj) += voigtnumber * ( StressTensor(k,l) * Identity(p,q) - StressTensor(l,q) * Identity( p, q)   );
               }
            }
         }
      }

      if ( rElementVariables.voigtsize == 6)
      {
         ThisMatrixSize = ThisMatrix;
      }
      else {
         for (unsigned int i = 0; i < rElementVariables.voigtsize; i++) {
            for (unsigned int j = 0; j < rElementVariables.voigtsize; j++) {
               unsigned int indexi = i; 
               unsigned int indexj = j;
               if ( i > 1) indexi += 1; 
               if ( j > 1) indexj += 1; 
               ThisMatrixSize(i, j) = ThisMatrix(indexi,indexj); 
            }
         }
      }

      MatrixType Kh=rLeftHandSideMatrix;

      Matrix Kuu = prod( trans( rVariables.B), rIntegrationWeight * Matrix( prod( ThisMatrixSize, rVariables.B) ) );
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
                  rLeftHandSideMatrix(indexi+2*i,indexj+2*j)+=Kuu(indexi,indexj);
                  indexj++;
               }
            }
            indexi++;
         }
      }

      if ( this->Id() == 0) {
       std::cout<<std::endl;
      std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;
      }


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      MatrixType Kh=rLeftHandSideMatrix;
      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int voigtsize = 3;
      if (dimension == 3) 
         voigtsize = 6;

      Vector IdentityVector = ZeroVector(voigtsize);
      for (unsigned int i = 0; i < dimension; i++)
         IdentityVector(i) = 1.0;


      Vector ThisVector = prod( trans( rVariables.B), IdentityVector);

      Matrix Kup = ZeroMatrix( number_of_nodes*dimension, number_of_nodes); 
      for (unsigned int i = 0; i < number_of_nodes*dimension; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            Kup(i,j) = ThisVector(i)*rVariables.N[j];
         }
      }
      Kup *= rIntegrationWeight; 

      unsigned int indexi = 0;
      unsigned int indexj = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int idim = 0; idim < dimension; idim ++)
         {
            indexj = 0;
            for (unsigned int j = 0; j < number_of_nodes; j++)
            {
               rLeftHandSideMatrix(indexi + 2*i, indexj + 3*j + 3 ) += Kup(indexi, indexj);
               indexj++;
            }
            indexi++;
         }

      }



      //std::cout<<std::endl;
      //std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;
      //std::cout << " KUP " << Kup << std::endl;

      KRATOS_CATCH( "" )
   }

   // *********************** KuJ TERMS ***********************************************
   // *********************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddKuJElemUJP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();

      Matrix Identity = ZeroMatrix(6, 1);
      for (unsigned int i = 0; i < dimension; i++)
         Identity(i, 0) = 1.0;

      ConstitutiveMatrix = prod(ConstitutiveMatrix, Identity);
      ConstitutiveMatrix *= rElementVariables.Alpha; 

      // Add (2alpha - 1) * SV
      for (unsigned int i = 0; i < 6; i++)
         ConstitutiveMatrix(i, 0) += ( 2.0* rElementVariables.Alpha - 1.0) * rVariables.StressVector(i);


      if ( rElementVariables.voigtsize < 6)
      {
         ConstitutiveMatrix(2, 0) = 0.0;
         ConstitutiveMatrix(4, 0) = 0.0;
         ConstitutiveMatrix(5, 0) = 0.0;
      }

      if ( this->Id() == 0) {
         std::cout << " THIS CONST MATRIX " << ConstitutiveMatrix << std::endl;
         std::cout << " IDENTITY " << Identity << std::endl;
      }

      // add the out of plane matrix 
      if ( rElementVariables.voigtsize < 6) { // ie. Plane Strain
         ConstitutiveMatrix(2, 0) -= rVariables.StressVector(2); 
         for (unsigned int i = 0; i < 2; i++)
            ConstitutiveMatrix(2,0) += rElementVariables.Alpha * rVariables.ConstitutiveMatrix(2, i) * Identity(i,0); 
      }

      ConstitutiveMatrix /= rElementVariables.NodalJacobian; 

      if ( this->Id() == 0 ) {
         std::cout << " AFTER OUT OF PLANE " << ConstitutiveMatrix << std::endl;
      }
      // MULTIPLY BY THE BETA DEVIATORIC
      Matrix BetaDeviatoric = ZeroMatrix(6,6);
      for (unsigned int i = 0; i < 6; i++)
         BetaDeviatoric(i,i) = 1.0;
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            BetaDeviatoric(i,j) -= rElementVariables.Beta; 
         }
      }

      ConstitutiveMatrix = prod( BetaDeviatoric, ConstitutiveMatrix);

      if (this->Id() == 0) {
         std::cout << " BETA DEVIATORIC " << BetaDeviatoric << std::endl;
         std::cout << " THETA BREVE DERIVATIVE " << ConstitutiveMatrix << std::endl;
      }

      Matrix SmallConstMatrix = ZeroMatrix( rElementVariables.voigtsize, 1);
      if ( rElementVariables.voigtsize == 6) {
         SmallConstMatrix = ConstitutiveMatrix;
      } 
      else {
         for (unsigned int i = 0; i < rElementVariables.voigtsize; i++) {
            unsigned int indexi = i;
            if ( indexi == 2) indexi += 1;
            SmallConstMatrix(i,0) = ConstitutiveMatrix(indexi, 0);
         }
      }

      if ( this->Id() == 0) {
         std::cout << " SMALL MATRIX " << SmallConstMatrix << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
      }

      Matrix  SmallMatrix = ZeroMatrix( dimension*number_of_nodes, number_of_nodes); //??

      SmallMatrix = prod( trans( rVariables.B), rIntegrationWeight * SmallConstMatrix );

      Matrix KuJ = ZeroMatrix( number_of_nodes*dimension, number_of_nodes);

      for (unsigned int i = 0; i < number_of_nodes*dimension; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            KuJ(i,j) += SmallMatrix( i,0) * rVariables.N[j];
         }
      }


      // ARA HE DE POSAR LA MATRIU AL SEU LLOC ( )
      MatrixType Kh=rLeftHandSideMatrix;
      unsigned int indexi = 0;
      unsigned int indexj = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int idim = 0; idim < dimension; idim ++)
         {
            indexj = 0;
            for (unsigned int j = 0; j < number_of_nodes; j++)
            {
               for (unsigned int jdim = 0; jdim< 1; jdim ++)
               {
                  rLeftHandSideMatrix(indexi + 2*i , indexj + 3*j + 2) += KuJ(indexi, indexj);
               }
            }
            indexi++;
         }

      }

      //std::cout<<std::endl;
      //std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )
   }

   // ******************** KJu term *******************************************************
   // *************************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddKJuElemUJP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      MatrixType Kh=rLeftHandSideMatrix;

      unsigned int indexJ = dimension;

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            int indexu = (dimension + 2) *j ;
            for (unsigned int dime = 0; dime < dimension; dime++)
            {  // AQUEST POT SER EL TERME QUE ESTÂ MALAMENT
               rLeftHandSideMatrix(indexJ, indexu + dime) -= rVariables.N[i] * rVariables.DN_DX( j, dime) * rIntegrationWeight ;
            }
         }
         indexJ += (dimension + 2);
      }


      //std::cout<<std::endl;
      //std::cout<<" KJu "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }


   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^KJJ term ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   void UpdatedLagrangianUJPElement::CalculateAndAddKJJElemUJP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData&  rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


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

            rLeftHandSideMatrix(indexpi,indexpj)  += consistent * rIntegrationWeight / rVariables.detH;
            indexpj += (dimension+2);
         }

         indexpi += (dimension + 2);
      }

      KRATOS_CATCH( "" )
   }


   // ****************************** KJp TERM ( is zero ) ********************************
   // ************************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddKJp (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      // LMV: No et liis
      // aquest és l'unic que pot ser zero,....  ( CASI CONVENÇUT)

      KRATOS_CATCH( "" )
   }



   // ******************************** TERM KpJ ****************************************
   // **********************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddKpJ (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)

   {
      // THIS TERM IS INCORRECT ( or SOMETHING)
      KRATOS_TRY

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      MatrixType Kh=rLeftHandSideMatrix;


      double KNumber = 0;

      for (unsigned int i = 0; i < dimension; i++) 
      {
         for (unsigned int j = 0; j< dimension; j++) {
            KNumber += ConstitutiveMatrix(i,j);
         }
      }
      KNumber *= rElementVariables.Alpha;

      for (unsigned int i = 0; i < dimension; i++)
      {
         KNumber += ( 2.0*rElementVariables.Alpha - 1.0) * rElementVariables.StressVectorEC(i);
      }

      if ( dimension == 2) {
         KNumber -= rVariables.StressVector(2);
         for (unsigned int i = 0; i < dimension; i++) {
            KNumber += rElementVariables.Alpha * ConstitutiveMatrix(2,i);
         }
      }

      KNumber *=  rElementVariables.Beta / rElementVariables.NodalJacobian; 
      KNumber *= rIntegrationWeight / rVariables.detH ;


      double consistent = 1.0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            consistent = 1.0/12.0;
            if ( i == j )
               consistent *= 2.0;
            rLeftHandSideMatrix( (dimension+2)*i +3, (dimension+2)*j + 2) -= KNumber * rVariables.N[i] * rVariables.N[j]; 
         }
      }

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY


      // per mi que està semi malament perquè em surt 0 en matlab i aquí em surt algun numero (1,3) diferent de zero, i no ho veig

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      MatrixType Kh=rLeftHandSideMatrix;

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix; 


      // 1. Definition of some tensors
      // 1.a Compress ConstitutiveMatrix
      Matrix ConstMatrix = ZeroMatrix(rElementVariables.voigtsize,rElementVariables.voigtsize);

      if ( rElementVariables.voigtsize == 6)
      {
         ConstMatrix = ConstitutiveMatrix;
      }
      else {
         int indexi, indexj; 
         for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
               indexi = i; indexj = j;
               if (indexi == 2)
                  indexi += 1;
               if (indexj == 2)
                  indexj += 1;
               ConstMatrix(i,j) = ConstitutiveMatrix(indexi, indexj) ;
            }
         }

      }

      // 1.b Definition of Identity and deviatoric 
      Matrix Identity = ZeroMatrix( 1, rElementVariables.voigtsize );
      for (unsigned int i = 0; i < dimension; i++) {
         Identity(0, i) = 1.0;
      }

      Matrix AlphaDeviatoricMatrix = ZeroMatrix(rElementVariables.voigtsize);
      for (unsigned int i = 0; i < dimension; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            AlphaDeviatoricMatrix(i,j) = -rElementVariables.Alpha;
         }
      }
      for (unsigned int i = 0; i < rElementVariables.voigtsize; i++) {
         AlphaDeviatoricMatrix(i,i) += 1.0;
      }


      // EVERTHING HERE.
      ConstMatrix = prod( Identity, ConstMatrix);
      ConstMatrix = prod( ConstMatrix, AlphaDeviatoricMatrix);


      //Geometric like terms;
      for (unsigned int i = 0; i < rElementVariables.voigtsize; i++) {
         ConstMatrix(0,i) += 2.0 * rElementVariables.StressVectorEC(i);
      }


      for (unsigned int i = 0; i < dimension; i++) // to compute the pseudo-pressure
      {
         for (unsigned int j = 0; j < dimension; j++) {
            ConstMatrix(0, j) -=  2.0 * rElementVariables.Alpha * rElementVariables.StressVectorEC(i);
         }
      }


      // OUT OF PLANE IN PS
      if ( dimension == 2) {
         for ( unsigned i = 0; i < rElementVariables.voigtsize; i++) {
            for (unsigned j = 0; j < rElementVariables.voigtsize; j++) {
               unsigned int indexi = i;
               if ( i > 1)
                  indexi += 1;
               ConstMatrix( 0,j) += rVariables.ConstitutiveMatrix(2, indexi) * AlphaDeviatoricMatrix(i,j);
            }
         }
      }


      if (this->Id()  == 0) {
         std::cout << " THIS DERIVATIVE:_ THIS CONST MATRIX. Kpu " << ConstMatrix * rElementVariables.Beta << std::endl;
         std::cout << " derivative of mean pressure respect displacement " << std::endl;
         //std::cout << " STRESS " << rVariables.StressVector << std::endl;
      }

      ConstMatrix *= rIntegrationWeight / rVariables.detH;
      ConstMatrix *= rElementVariables.Beta; 
      ConstMatrix = prod( ConstMatrix, rVariables.B);

      // Multiply and everithing

      if ( this->Id() == 0 )
      {
         std::cout << " CONSMAT  SIZE    " << ConstMatrix.size1() << " and " << ConstMatrix.size2() << std::endl;
         std::cout << " B SIZE " << rVariables.B.size1() << " and " << rVariables.B.size2() << std::endl;
         std::cout << " JUST TO CHECK " << ConstMatrix.size2() << " and is equal to " << dimension*number_of_nodes << std::endl;
      }


      Matrix SmallMatrix = ZeroMatrix( number_of_nodes, dimension*number_of_nodes);
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes * dimension; j++) {
            SmallMatrix(i,j) += rVariables.N[i]* ConstMatrix(0,j); // NOT SURE
         }
      }

      unsigned int indexj = 0;
      for (unsigned int i = 0; i < number_of_nodes ; i++)
      {
         indexj = 0;
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            for (unsigned jdim = 0; jdim<dimension; jdim++) {
               rLeftHandSideMatrix( dimension + 1 + i*(2+dimension), j*(2+dimension) + jdim ) -= SmallMatrix(i,indexj);
               if ( this->Id() == 0)
                  std::cout << "CHECK INDICES: TO " << dimension + 1 + i*(2+dimension) << " and " << (2+dimension)*j + jdim << "FROM " << i << " and " << indexj << std::endl;
               indexj += 1;
            }
         }
      }

      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension+1;

      double consistent ;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension+1;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent = 1.0/12.0;
            if ( i == j)
               consistent *= 2.0;
            rLeftHandSideMatrix(indexpi,indexpj)  += consistent * rIntegrationWeight / rVariables.detH;
            indexpj += (dimension + 2);
         }

         indexpi += (dimension + 2);
      }

      KRATOS_CATCH( "" )
   }



   //************************************************************************************
   //************************************************************************************
   void UpdatedLagrangianUJPElement::CalculateAndAddKJJStabElemUJP (MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         UJPElementData & rElementVariables, 
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

      double AlphaStabilization  = 4.0; 
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
      AlphaStabilization *= StabilizationFactor; 

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

      AlphaStabilization *= BulkModulus;  // TIMES THE BULK MODULUS BECAUSE I HAVE ALL THE EQUATION MULTIPLIED BY THE BULK MODULUS
      if ( YoungModulus < 0.0001)
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0;
         AlphaStabilization *= mElementStabilizationNumber; 

      }

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension ;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=(-1)*AlphaStabilization;
            if(indexpi==indexpj)
               consistent=2*AlphaStabilization;

            rLeftHandSideMatrix(indexpi,indexpj) += consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);     //2D

            indexpj += (dimension + 2);
         }

         indexpi += (dimension + 2);
      } 
      KRATOS_CATCH( "" )
   }

   void UpdatedLagrangianUJPElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         UJPElementData & rElementVariables, 
         double& rIntegrationWeight)
   {

      KRATOS_TRY

      //repasar

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension+1;
      double consistent = 1.0;

      double AlphaStabilization  = 4.0; 
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_P];
      AlphaStabilization *= StabilizationFactor;

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2.0*(1.0+PoissonCoefficient));
      double Bulk = YoungModulus/ ( 3.0 * ( 1.0 - 2.0*PoissonCoefficient) );

      AlphaStabilization = (4.0*StabilizationFactor/(18.0*LameMu)) * Bulk;  // TIMES THE BULK MODULUS
      if ( YoungModulus < 0.0001)
         return;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension + 1;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=(-1.0)*AlphaStabilization;
            if(indexpi==indexpj)
               consistent=2.0*AlphaStabilization;

            rLeftHandSideMatrix(indexpi,indexpj) += consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);     //2D

            indexpj += (dimension + 2);
         }

         indexpi += (dimension + 2);
      }

      KRATOS_CATCH( "" )
   }





   // ^****************** CalculateThisElementVariables ******************************************
   // *********** Compute only once some terms **********************************************************
   void UpdatedLagrangianUJPElement::CalculateThisElementData( UJPElementData& rElementVariables, const ElementDataType & rVariables)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      rElementVariables.voigtsize = 3;
      rElementVariables.Beta = 1.0/3.0;
      rElementVariables.Alpha = 0.5;
      if ( dimension == 3) {
         rElementVariables.voigtsize = 6;
         rElementVariables.Alpha = 1.0/3.0;
      }


      rElementVariables.NodalMeanStress = 0;
      rElementVariables.NodalJacobian = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         rElementVariables.NodalMeanStress += GetGeometry()[i].GetSolutionStepValue( PRESSURE ) * rVariables.N[i];
         rElementVariables.NodalJacobian   += GetGeometry()[i].GetSolutionStepValue( JACOBIAN ) * rVariables.N[i];
      }

      rElementVariables.ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3 ; i++)
         rElementVariables.ElementalMeanStress += rVariables.StressVector(i);
      rElementVariables.ElementalMeanStress /= 3.0;

      Vector AuxStress = ZeroVector(6);
      AuxStress = rVariables.StressVector; 
      for (unsigned int i = 0; i < 3; i++)
         AuxStress(i) += ( rElementVariables.NodalMeanStress - rElementVariables.ElementalMeanStress);

      rElementVariables.StressVector = ZeroVector(rElementVariables.voigtsize);
      rElementVariables.StressVectorEC = ZeroVector( rElementVariables.voigtsize );
      if ( rElementVariables.voigtsize == 6) {
         rElementVariables.StressVector = AuxStress; 
         rElementVariables.StressVectorEC = rVariables.StressVector; 
      }
      else {
         rElementVariables.StressVector(0) = AuxStress(0);
         rElementVariables.StressVector(1) = AuxStress(1);
         rElementVariables.StressVector(2) = AuxStress(3);

         rElementVariables.StressVectorEC(0) = rVariables.StressVector(0);
         rElementVariables.StressVectorEC(1) = rVariables.StressVector(1);
         rElementVariables.StressVectorEC(2) = rVariables.StressVector(3);
      }

   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void UpdatedLagrangianUJPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }



}


