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
#include "custom_elements/updated_lagrangian_U_P_Second_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

// **
// UP ELEMENT BUT NOT FOR ISOCHORIC PLASTICITY THAT HAS AN HYPERELASTIC SPLIT BETWEEN VOLUMETRIC AND DEVIATORIC PART
// ( i.e: it is for the MCC, but works for everything)
// VALE. CREC Q TINC LES TANGENTS BÉ. (peró no se sap mai)
// **

// OBSERVATION: I compute the constitutive equation with the appropiate law ( PlaneStrain, ...) but I get the 3D stress tensor or 3D constitutive matrix.
// This is why I have to "reform" the vector of matrix here.
// ( I do that compute the appropiate 3D pressure and not the pressure in the plane)

// he de parar d'intentar programar això i fer unes GeneralVariables per no haver de calcular coses varies vegades.

// VALEEEEEEE: Crec que ja tinc tot el tema de la Voigt notation bé, em sembla, però bueno, tot és possible


// Pot encara tenir un error, pq amb les lleis de JMC en versió 3D no acaba de convergir del tot i en canvi la UJP si que ho fa

// VALE, hi ha la mElementScalingNumber, que, de moment, hi poso 1/K perquè llavors convergeix millor perquè multiplico la segona equació per això.

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUPSecondElement::UpdatedLagrangianUPSecondElement()
      : UpdatedLagrangianUPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPSecondElement::UpdatedLagrangianUPSecondElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPSecondElement::UpdatedLagrangianUPSecondElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUPSecondElement::UpdatedLagrangianUPSecondElement( UpdatedLagrangianUPSecondElement const& rOther)
      :UpdatedLagrangianUPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUPSecondElement&  UpdatedLagrangianUPSecondElement::operator=(UpdatedLagrangianUPSecondElement const& rOther)
   {
      UpdatedLagrangianUPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUPSecondElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUPSecondElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUPSecondElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUPSecondElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUPSecondElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUPSecondElement::~UpdatedLagrangianUPSecondElement()
   {
   }



   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUPSecondElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);

      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P element type ", " UpdatedLagrangianUPElement" );

      //verify that the variables are correctly initialized

      if ( PRESSURE.Key() == 0 )
         KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" );


      return correct;

      KRATOS_CATCH( "" );
   }

   //*************************************************************************
   //*************************************************************************

   void UpdatedLagrangianUPSecondElement::GetValueOnIntegrationPoints( const Variable<Matrix >& rVariable, 
         std::vector<Matrix>& rValues,
         const ProcessInfo& rCurrentProcessInfo)
   {

      KRATOS_TRY
      if ( rVariable == EQ_CAUCHY_STRESS) {
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

   void UpdatedLagrangianUPSecondElement::CalculateOnIntegrationPoints( const Variable <Matrix >& rVariable, std::vector<Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rOutput.size() != integration_points_number )
         rOutput.resize( integration_points_number );

      if ( rVariable == EQ_CAUCHY_STRESS)
      {
         const unsigned int number_of_nodes = GetGeometry().PointsNumber();

         //create and initialize element variables:
         GeneralVariables Variables;
         this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

         //create constitutive law parameters:
         ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

         //set constitutive law flags:
         Flags &ConstitutiveLawOptions = Values.GetOptions();

         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

         // for integration points
         for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
         {

            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //to take in account previous step writing
            if( mFinalizedStep ){
               this->GetHistoricalVariables(Variables,PointNumber);
            }		

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

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
      else {
         LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
         //UpdatedLagrangianUPElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }


      KRATOS_CATCH( "" )
   }

   //********************** COMPUTE LHS **********************************************************
   //********************** COMPUTE LHS **********************************************************
   void UpdatedLagrangianUPSecondElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
   {

      // Scaling constant for the second equation.
      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      mElementScalingNumber = 1.0;
      if (YoungModulus > 1e-5 )
         mElementScalingNumber = 100.0/BulkModulus; 

      
      
      
      //contributions to stiffness matrix calculated on the reference config
      rVariables.detF0   *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // operation performed: add Km to the rLefsHandSideMatrix

      // simplifies the code and some things are only computed once.
      ThisElementGeneralVariables ElementVariables;
      CalculateThisElementGeneralVariables( ElementVariables, rVariables);

      //respect to the current configuration n+1
      CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kg to the rLefsHandSideMatrix
      CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kup to the rLefsHandSideMatrix
      CalculateAndAddKup( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpu to the rLefsHandSideMatrix
      CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpp to the rLefsHandSideMatrix
      CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      // operation performed: add Kpp Stab to the rLefsHandSideMatrix
      CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );

      //KRATOS_WATCH( rLeftHandSideMatrix )

      rVariables.detF     = DeterminantF;
      rVariables.detF0   /= rVariables.detF;
      //KRATOS_WATCH( rLeftHandSideMatrix )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPSecondElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      // Scaling constant for the second equation.
      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      mElementScalingNumber = 1.0;
      if (YoungModulus > 1e-5 )
         mElementScalingNumber = 100.0/BulkModulus; 

      
      
      //contribution to external forces
      rVariables.detF0   *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

      // simplifies the code and some things are only computed once.
      ThisElementGeneralVariables ElementVariables;
      CalculateThisElementGeneralVariables( ElementVariables, rVariables);

      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
      CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

      // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
      CalculateAndAddInternalForces( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
      CalculateAndAddPressureForces( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

      //KRATOS_WATCH( rRightHandSideVector )

      rVariables.detF     = DeterminantF;
      rVariables.detF0   /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }

   //*************************************************************************
   //*************************************************************************
   void UpdatedLagrangianUPSecondElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
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
         unsigned int indexup = dimension * i + i;
         unsigned int indexu  = dimension * i;

         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
         }
      }

      /*if (this->Id() < 10) {
        std::cout << " THIS IS THE STRESS in InternalForces " << rVariables.StressVector << std::endl;
        std::cout << " FF " << rVariables.F << std::endl;
        std::cout << " FF2 " << rVariables.FT << std::endl;
        }*/

      KRATOS_CATCH( "" )
   }
   //*************************************************************************
   //*************************************************************************

   void UpdatedLagrangianUPSecondElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;
      VectorType Fh=rRightHandSideVector;
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

   void UpdatedLagrangianUPSecondElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
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
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
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
         //this->GetValueOnIntegrationPoints( SIMILAR_SHEAR_MODULUS, Values, SomeProcessInfo);
         LargeDisplacementElement::GetValueOnIntegrationPoints( SIMILAR_SHEAR_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization /= Values[0];

         //GetValueOnIntegrationPoints( SIMILAR_BULK_MODULUS, Values, SomeProcessInfo);
         LargeDisplacementElement::GetValueOnIntegrationPoints( SIMILAR_BULK_MODULUS, Values, SomeProcessInfo);
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

            // std::cout<<" Pressure "<<Pressure<<std::endl;
         }


         indexp += (dimension + 1);
      }



      KRATOS_CATCH( "" )

   }

   // ******************************** KUUM **************************************************
   // ****************************************************************************************
   void UpdatedLagrangianUPSecondElement::CalculateAndAddKuum ( MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      Matrix ConstitutiveMatrix = ZeroMatrix(rElementVariables.voigtsize); 

      ECConstitutiveMatrix = prod( rElementVariables.DeviatoricTensor, ECConstitutiveMatrix);

      // that collection of terms.

      if ( this->Id() < 0) {
         std::cout << " CONST MATRIX FROM EECC " << ECConstitutiveMatrix << std::endl;
      }


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
           ECConstitutiveMatrix(i,j) -= (2.0/3.0) *voigtnumber* rVariables.StressVector(j);
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
   void UpdatedLagrangianUPSecondElement::CalculateAndAddKuug ( MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      int size = number_of_nodes * dimension;

      Vector StressVector = rElementVariables.StressVector;

      Matrix StressTensor = MathUtils<double>::StressVectorToTensor( StressVector );

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

   void UpdatedLagrangianUPSecondElement::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         ThisElementGeneralVariables& rElementVariables, 
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

      // std::cout<<std::endl;
      // std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPSecondElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         ThisElementGeneralVariables& rElementVariables, 
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

      if ( this->Id() < 0)
      {
         std::cout << " THIS DERIVATIVE IS " << ConstitutiveVector << std::endl; 
         double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
         std::cout << " JMC " << BulkModulus * ( 1 - log(rVariables.detFT) ) / rVariables.detFT << std::endl;
         std::cout << " J " << rVariables.detFT << std::endl;
            std::cout << std::endl;

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

   void UpdatedLagrangianUPSecondElement::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         ThisElementGeneralVariables& rElementVariables, 
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
   void UpdatedLagrangianUPSecondElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
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
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
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
         //this->GetValueOnIntegrationPoints( SIMILAR_SHEAR_MODULUS, Values, SomeProcessInfo);
         LargeDisplacementElement::GetValueOnIntegrationPoints( SIMILAR_SHEAR_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization /= Values[0];

         //GetValueOnIntegrationPoints( SIMILAR_BULK_MODULUS, Values, SomeProcessInfo);
         LargeDisplacementElement::GetValueOnIntegrationPoints( SIMILAR_BULK_MODULUS, Values, SomeProcessInfo);
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



   // CHARACTERISTIC ELEMENT SIZE ACCORDING TO SOME PAPER
   double UpdatedLagrangianUPSecondElement::GetElementSize( const Matrix& rDN_DX)
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
      // this size is recomended from Sun, Ostien and Salinger (IJNAMG, 2013)
      // the same multiplied by two in order to...
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPSecondElement::SetGeneralVariables(GeneralVariables& rVariables,
         ConstitutiveLaw::Parameters& rValues,
         const int & rPointNumber)
   {

      if(rVariables.detF<0){

         std::cout<<" Element: "<<this->Id()<<std::endl;

         unsigned int number_of_nodes = GetGeometry().PointsNumber();

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
            std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
         }

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
               array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
               array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
               std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
            }
            else{
               std::cout<<" ---Contact_Force: NULL "<<std::endl;
            }
         }
         std::cout << " detF " << rVariables.detF << " rVariables.detF0 " << rVariables.detF0 << std::endl;
         std::cout << " F " << rVariables.F << " F0 " << rVariables.F0 << std::endl;
         KRATOS_THROW_ERROR( std::invalid_argument," LARGE DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
      }

      //Compute total F: FT 
      rVariables.detFT = rVariables.detF * rVariables.detF0;
      rVariables.FT    = prod( rVariables.F, rVariables.F0 );

      rValues.SetDeterminantF(rVariables.detFT);
      rValues.SetDeformationGradientF(rVariables.FT);
      rValues.SetStrainVector(rVariables.StrainVector);
      rValues.SetStressVector(rVariables.StressVector);
      rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
      rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
      rValues.SetShapeFunctionsValues(rVariables.N);
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPSecondElement::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
   {

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      unsigned int voigtsize  = 3;

      if( dimension == 3 )
      {
         voigtsize  = 6;
      }

      rVariables.detF  = 1;

      rVariables.detF0 = 1;

      rVariables.detFT = 1;

      rVariables.detJ = 1;

      rVariables.B.resize( voigtsize , number_of_nodes * dimension );

      rVariables.F.resize( dimension, dimension );

      rVariables.F0.resize( dimension, dimension );

      rVariables.FT.resize( dimension, dimension );

      rVariables.ConstitutiveMatrix.resize( 6, 6); //  (voigtsize, voigtsize );

      rVariables.StrainVector.resize( voigtsize );

      rVariables.StressVector.resize( 6 ); // try to see what happens...

      rVariables.DN_DX.resize( number_of_nodes, dimension );

      //set variables including all integration points values

      //reading shape functions
      rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

      //reading shape functions local gradients
      rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

      //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
      rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

      // AND SOME PART THAT I MISSED
      //Calculate Delta Position
      rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

      //set variables including all integration points values

      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );



   }

   void UpdatedLagrangianUPSecondElement::CalculateThisElementGeneralVariables( ThisElementGeneralVariables& rElementGeneralVariables, const GeneralVariables & rVariables)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      rElementGeneralVariables.voigtsize = 3;
      if ( dimension == 3)
         rElementGeneralVariables.voigtsize = 6;

      rElementGeneralVariables.NodalMeanStress = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         rElementGeneralVariables.NodalMeanStress += GetGeometry()[i].GetSolutionStepValue( PRESSURE ) * rVariables.N[i];

      rElementGeneralVariables.ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3 ; i++)
         rElementGeneralVariables.ElementalMeanStress += rVariables.StressVector(i);
      rElementGeneralVariables.ElementalMeanStress /= 3.0;

      Vector AuxStress = ZeroVector(6);
      AuxStress = rVariables.StressVector; 
      for (unsigned int i = 0; i < 3; i++)
         AuxStress(i) += ( rElementGeneralVariables.NodalMeanStress - rElementGeneralVariables.ElementalMeanStress);

      rElementGeneralVariables.StressVector = ZeroVector(rElementGeneralVariables.voigtsize);

      if ( rElementGeneralVariables.voigtsize == 6) {
         rElementGeneralVariables.StressVector = AuxStress; 
      }
      else {
         rElementGeneralVariables.StressVector(0) = AuxStress(0);
         rElementGeneralVariables.StressVector(1) = AuxStress(1);
         rElementGeneralVariables.StressVector(2) = AuxStress(3);
      }

      rElementGeneralVariables.DeviatoricTensor = ZeroMatrix(6);
      for (unsigned int i = 0; i < 6; i++)
         rElementGeneralVariables.DeviatoricTensor(i,i) = 1;

      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
         rElementGeneralVariables.DeviatoricTensor(i,j) -= 1.0/3.0;
         }
      }



   }
      

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUPSecondElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void UpdatedLagrangianUPSecondElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }



}  // Namespace Kratos
