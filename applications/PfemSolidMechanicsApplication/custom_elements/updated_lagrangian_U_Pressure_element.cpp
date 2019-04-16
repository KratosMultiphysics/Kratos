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
#include "custom_elements/updated_lagrangian_U_Pressure_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"

// **
// UP ELEMENT BUT NOT FOR ISOCHORIC PLASTICITY THAT HAS AN HYPERELASTIC SPLIT BETWEEN VOLUMETRIC AND DEVIATORIC PART
// ( i.e: it is for the MCC, but works for everything)
// **

// OBSERVATION: I compute the constitutive equation with the appropiate law ( PlaneStrain, ...) but I get the 3D stress tensor or 3D constitutive matrix.
// This is why I have to "reform" the vector of matrix here.


namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUPressureElement::UpdatedLagrangianUPressureElement()
      : UpdatedLagrangianUPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPressureElement::UpdatedLagrangianUPressureElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPressureElement::UpdatedLagrangianUPressureElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUPressureElement::UpdatedLagrangianUPressureElement( UpdatedLagrangianUPressureElement const& rOther)
      :UpdatedLagrangianUPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUPressureElement&  UpdatedLagrangianUPressureElement::operator=(UpdatedLagrangianUPressureElement const& rOther)
   {
      UpdatedLagrangianUPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUPressureElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUPressureElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUPressureElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUPressureElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUPressureElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPressureElement::~UpdatedLagrangianUPressureElement()
   {
   }



   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUPressureElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      //correct = LargeDisplacementElement::Check(rCurrentProcessInfo);
      // check abreviat. Intento fer servir una llei 3D en tots els casos, a veure què passa


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P element type ", " UpdatedLagrangianUPElement" );

      //verify that the variables are correctly initialized

      if ( PRESSURE.Key() == 0 )
         KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" );

      if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW)->GetStrainSize() != 6) {
         KRATOS_THROW_ERROR( std::invalid_argument, " since I do not know how to do it correctly, I try to have a 3D law in here ", this->Id() );
      }

      unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
      if ( dimension == 2 )
      {
         if ( this->GetProperties().Has( THICKNESS ) == false ){
            this->GetProperties().SetValue( THICKNESS , 1.0 );
         }
      }


      return correct;

      KRATOS_CATCH( "" );
   }

   //*************************************************************************
   //*************************************************************************

   void UpdatedLagrangianUPressureElement::GetValueOnIntegrationPoints( const Variable<Matrix >& rVariable, 
         std::vector<Matrix>& rValues,
         const ProcessInfo& rCurrentProcessInfo)
   {

      KRATOS_TRY
      if ( rVariable == CAUCHY_STRESS_TENSOR) {
         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);
      }
      else if ( rVariable == TOTAL_CAUCHY_STRESS) {
         CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);
      }
      else {
         LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);
         //UpdatedLagrangianUPElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);
      }

      KRATOS_CATCH( "" )

   }

   //*************************************************************************
   //*************************************************************************

   void UpdatedLagrangianUPressureElement::GetValueOnIntegrationPoints( const Variable<double> & rVariable,
         std::vector<double>& rValues,
         const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      //Element* pBasePointer = this;
      //pBasePointer->GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);
      SolidElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);

      KRATOS_CATCH("")
   }

   //*************************************************************************
   //*************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateOnIntegrationPoints( const Variable <Matrix >& rVariable, std::vector<Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rOutput.size() != integration_points_number )
         rOutput.resize( integration_points_number );

      if ( rVariable == CAUCHY_STRESS_TENSOR)
      {
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

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);

            Vector StressVector = Variables.StressVector;
            double ElementalPressure = 0;
            double NodalPressure = 0;
            for (unsigned int i = 0; i < number_of_nodes; i++)
               NodalPressure += GetGeometry()[i].GetSolutionStepValue(PRESSURE) * Variables.N[i];
            for (unsigned int i = 0; i < 3; i++)
               ElementalPressure += StressVector(i);
            ElementalPressure /= 3.0; 

            for (unsigned int i = 0; i < 3; i++)
               StressVector(i) += ( NodalPressure - ElementalPressure);


            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector); 

         }

      }
      else if ( rVariable == TOTAL_CAUCHY_STRESS) {

         CalculateOnIntegrationPoints( CAUCHY_STRESS_TENSOR, rOutput, rCurrentProcessInfo);

         if ( GetGeometry()[0].HasDofFor( WATER_PRESSURE) ) 
         {
            const unsigned int number_of_nodes = GetGeometry().size();
            //create and initialize element variables:
            ElementDataType Variables;
            this->InitializeElementData(Variables,rCurrentProcessInfo);

            //reading integration points
            for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
            {
               //compute element kinematics B, F, DN_DX ...
               this->CalculateKinematics(Variables,PointNumber);

               double WaterPressure = 0;
               for (unsigned int i = 0; i < number_of_nodes; i++)
               {
                  WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE) * Variables.N[i];
               }

               for (unsigned int i = 0; i < 3; i++)
                  rOutput[PointNumber](i,i) += WaterPressure; 

            }
         }
      }
      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
         //UpdatedLagrangianUPElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }


      KRATOS_CATCH( "" )
   }

   //********************** COMPUTE LHS **********************************************************
   //********************************************************************************
   void UpdatedLagrangianUPressureElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      KRATOS_TRY

      // Scaling constant for the second equation.
      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      mElementScalingNumber = 1.0;

      if (YoungModulus > 1e-5 ) {
         mElementScalingNumber = 100.0/BulkModulus; 
      }
      else {
         mElementScalingNumber = 0.01;
      }
      mElementScalingNumber = 1.0;
      
      //contributions to stiffness matrix calculated on the reference config
      rVariables.detF0   *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // operation performed: add Km to the rLefsHandSideMatrix

      // simplifies the code and some things are only computed once.
      ThisElementData ElementVariables;
      this->CalculateThisElementData( ElementVariables, rVariables);

      //respect to the current configuration n+1
      this->CalculateAndAddKuumElemUP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kg to the rLefsHandSideMatrix
      this->CalculateAndAddKuugElemUP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kup to the rLefsHandSideMatrix
      this->CalculateAndAddKupElemUP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpu to the rLefsHandSideMatrix
      this->CalculateAndAddKpuElemUP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpp to the rLefsHandSideMatrix
      this->CalculateAndAddKppElemUP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpp Stab to the rLefsHandSideMatrix
      this->CalculateAndAddKppStabElemUP( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );


      rVariables.detF     = DeterminantF;
      rVariables.detF0   /= rVariables.detF;

      KRATOS_CATCH("")
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      KRATOS_TRY

      // Scaling constant for the second equation.
      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      mElementScalingNumber = 1.0;

      if (YoungModulus > 1e-5 ) {
         mElementScalingNumber = 100.0/BulkModulus; 
      }
      else {
         mElementScalingNumber = 0.01;
      }
      mElementScalingNumber = 1.0;
      
      
      //contribution to external forces
      rVariables.detF0   *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

      // simplifies the code and some things are only computed once.
      ThisElementData ElementVariables;
      CalculateThisElementData( ElementVariables, rVariables);

      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
      CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

      // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
      CalculateAndAddInternalForcesElemUP( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
      CalculateAndAddPressureForcesElemUP( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedPressureElemUP( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      rVariables.detF     = DeterminantF;
      rVariables.detF0   /= rVariables.detF;

      KRATOS_CATCH("")
   }

   //*************************************************************************
   //*************************************************************************
   void UpdatedLagrangianUPressureElement::CalculateAndAddInternalForcesElemUP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //VectorType Fh=rRightHandSideVector;

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rElementVariables.StressVector );

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
   //*************************************************************************
   //*************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateAndAddPressureForcesElemUP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;
      //VectorType Fh=rRightHandSideVector;
      double consistent=1;

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

         rRightHandSideVector[indexp] -= rElementVariables.ElementalMeanStress * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;


         indexp += (dimension + 1);
      }


      KRATOS_CATCH( "" )

   }



   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateAndAddStabilizedPressureElemUP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

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


         indexp += (dimension + 1);
      }



      KRATOS_CATCH( "" )

   }

   // ******************************** KUUM **************************************************
   // ****************************************************************************************
   void UpdatedLagrangianUPressureElement::CalculateAndAddKuumElemUP ( MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      Matrix ConstitutiveMatrix = ZeroMatrix(rElementVariables.voigtsize,rElementVariables.voigtsize); 

      ECConstitutiveMatrix = prod( rElementVariables.DeviatoricTensor, ECConstitutiveMatrix);

      // that collection of terms.

      // a. 2 ( pElem - p Nodal) I4S
      for ( unsigned int i = 0; i < 6; i++) {
         double voigtnumber = 1.0;
         if ( i > 2)
            voigtnumber = 0.500;
         ECConstitutiveMatrix(i,i) += 2.0*voigtnumber * ( rElementVariables.ElementalMeanStress - rElementVariables.NodalMeanStress);
      }

      // b. p 1 cross 1
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            ECConstitutiveMatrix(i,j) += rElementVariables.NodalMeanStress ;
         }
      }

      // c . - 2/3 I times sigma
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 6; j++) {
           double voigtnumber = 1.0;
           if ( j > 2)
              voigtnumber = 1.00;
           ECConstitutiveMatrix(i,j) -= (2.0/6.0) *voigtnumber* rVariables.StressVector(j);
         }
      }

      // LMV: try to put the sim part
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 6; j++) {
           double voigtnumber = 1.0;
           if ( j > 2)
              voigtnumber = 1.00;
           ECConstitutiveMatrix(j,i) -= (2.0/6.0) *voigtnumber* rVariables.StressVector(j);
         }
      }


      if ( rElementVariables.voigtsize == 6)
      {
         ConstitutiveMatrix = ECConstitutiveMatrix;
      }
      else
      {
         int indexi , indexj; 
         for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3 ; j++)
            {
               indexi = i;
               indexj = j;
               if ( i == 2)
                  indexi += 1;
               if ( j == 2)
                  indexj += 1;

               ConstitutiveMatrix(i,j) = ECConstitutiveMatrix(indexi, indexj);
            }
         }
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
   void UpdatedLagrangianUPressureElement::CalculateAndAddKuugElemUP ( MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      int size = number_of_nodes * dimension;

      Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rElementVariables.StressVector );

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

   //********************* KuP ***********************************************
   //*************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateAndAddKupElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         ThisElementData& rElementVariables, 
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
               rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight; // * rVariables.detF; (aquest és el que posem a u, o sigui que no cal comentar)
            }
            indexp += (dimension + 1);
         }
      }

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateAndAddKpuElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      //repasar

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      Matrix Kpu;
      //MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration


      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;

      Matrix Identity = ZeroMatrix(1,rElementVariables.voigtsize);
      for (unsigned int i = 0; i < 3 ; ++i) {
         Identity(0,i) = 1.0; 
      }

      Matrix AuxConstitutiveVector = prod( Identity, ECConstitutiveMatrix);
      AuxConstitutiveVector *= ( 1.0/ 3.0 );

      // EXTRA TERMS.
      // 2 / 3 * Stress
      for (unsigned int i = 0; i < 6; i++) {
         double voigtnumber = 1.0;
         if ( i > 2)
            voigtnumber = 1.0; // ithink,...
         AuxConstitutiveVector(0, i) += (2.0/3.0) * voigtnumber * rVariables.StressVector(i);
      }
      // - ElemenPressure * 1
      for (unsigned int i = 0; i < 3; i++)
         AuxConstitutiveVector(0,i) -= rElementVariables.ElementalMeanStress;


      Matrix ConstitutiveVector = ZeroMatrix( 1, rElementVariables.voigtsize);
      if ( rElementVariables.voigtsize == 6)
      {
         ConstitutiveVector = AuxConstitutiveVector;
      }
      else
      {
         int indexi; 
         for (unsigned int i = 0; i < 3; i++) {
            indexi = i;
            if ( i == 2)
               indexi += 1;

            ConstitutiveVector(0,i) = AuxConstitutiveVector(0,indexi);
         }
      }

        Kpu = prod( ConstitutiveVector, rVariables.B);

      Matrix Kpu2 = ZeroMatrix( number_of_nodes, number_of_nodes*dimension);
      for (unsigned int i = 0; i < number_of_nodes;  i++) {
         for (unsigned int j = 0; j < number_of_nodes*dimension; j++) {
            Kpu2 ( i, j ) = rVariables.N[i] * Kpu(0,j);
         }
      }
      Kpu2 *= rIntegrationWeight;

      unsigned int indexp = dimension;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexp,indexup+k) +=  Kpu2( i, j*2 + k) *mElementScalingNumber;
            }
         }
         indexp += (dimension + 1);
      }


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPressureElement::CalculateAndAddKppElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         ThisElementData& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;
      double consistent = 1.0;


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=1;
            if(indexpi==indexpj)
               consistent=2;

            rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * (1.0/12.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber; //2D

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      // std::cout<<std::endl;
      // std::cout<<" Kpp "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }





   //************************************************************************************
   //************************************************************************************
   void UpdatedLagrangianUPressureElement::CalculateAndAddKppStabElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables, 
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
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_P];
      AlphaStabilization *= StabilizationFactor;
      if ( StabilizationFactor < 0.0001)
         return;

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      //Experimental
      // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
      //   LameMu = rVariables.ConstitutiveMatrix(2,2);


      //use of this variable for the complete parameter:
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

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=(-1)*AlphaStabilization;
            if(indexpi==indexpj)
               consistent=2*AlphaStabilization;

            rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;     //2D

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

   void UpdatedLagrangianUPressureElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
   {

      UpdatedLagrangianUPElement::InitializeElementData(rVariables, rCurrentProcessInfo);

      // calculate the constitutive equation as a 3D.
      rVariables.ConstitutiveMatrix.resize( 6, 6); //  (voigtsize, voigtsize );

      rVariables.StressVector.resize( 6 ); // try to see what happens...

   }

   void UpdatedLagrangianUPressureElement::CalculateThisElementData( ThisElementData& rElementVariables, const ElementDataType & rVariables)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      rElementVariables.voigtsize = 3;
      if ( dimension == 3)
         rElementVariables.voigtsize = 6;

      rElementVariables.NodalMeanStress = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         rElementVariables.NodalMeanStress += GetGeometry()[i].GetSolutionStepValue( PRESSURE ) * rVariables.N[i];

      rElementVariables.ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3 ; i++)
         rElementVariables.ElementalMeanStress += rVariables.StressVector(i);
      rElementVariables.ElementalMeanStress /= 3.0;

      Vector AuxStress = ZeroVector(6);
      AuxStress = rVariables.StressVector; 
      for (unsigned int i = 0; i < 3; i++)
         AuxStress(i) += ( rElementVariables.NodalMeanStress - rElementVariables.ElementalMeanStress);

      rElementVariables.StressVector = ZeroVector(rElementVariables.voigtsize);

      if ( rElementVariables.voigtsize == 6) {
         rElementVariables.StressVector = AuxStress; 
      }
      else {
         rElementVariables.StressVector(0) = AuxStress(0);
         rElementVariables.StressVector(1) = AuxStress(1);
         rElementVariables.StressVector(2) = AuxStress(3);
      }

      rElementVariables.DeviatoricTensor = ZeroMatrix(6,6);
      for (unsigned int i = 0; i < 6; i++)
         rElementVariables.DeviatoricTensor(i,i) = 1;

      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
         rElementVariables.DeviatoricTensor(i,j) -= 1.0/3.0;
         }
      }



   }
      

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPressureElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUPElement )
   }

   void UpdatedLagrangianUPressureElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUPElement )
   }



}  // Namespace Kratos
