//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:          October 2017 $
//   Revision:            $Revision:               0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/sanisand_small_strain_model.hpp"

/* WRAPPER */
extern "C" void umat_wrapper_( double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
      double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
      double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
      char* MATERL, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
      double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
      double** DFGRD1, int* NOEL, int* NPT, double* KSLAY, double* KSPT, double* KSTEP,
      double* KINC, int* MATERIALNUMBER );


namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SanisandSmallStrainUmatModel::SanisandSmallStrainUmatModel()
      : SmallStrainUmatModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SanisandSmallStrainUmatModel::SanisandSmallStrainUmatModel(const SanisandSmallStrainUmatModel& rOther)
      : SmallStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer SanisandSmallStrainUmatModel::Clone() const
   {
      return Kratos::make_shared<SanisandSmallStrainUmatModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   SanisandSmallStrainUmatModel& SanisandSmallStrainUmatModel::operator=(SanisandSmallStrainUmatModel const& rOther)
   {
      SmallStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SanisandSmallStrainUmatModel::~SanisandSmallStrainUmatModel()
   {
   }

   //************************************************************************************
   //************************************************************************************

   void SanisandSmallStrainUmatModel::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
   {
      KRATOS_TRY


      UmatDataType Variables;
      this->InitializeElasticData(rValues,Variables);

      // allocate variables
      int ndi = 3, nshr = 3, ntens = 6; // number of stress and strain components

      double pConstitutiveMatrix[ntens][ntens];
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            pConstitutiveMatrix[i][j] = 0;
         }
      }

      // time
      const ModelDataType& rModelData = Variables.GetModelData();
      double pTime[2];
      double delta_time = rModelData.GetProcessInfo()[DELTA_TIME];
      if ( delta_time < 1e-5)
         delta_time = 0.01;
      pTime[0] = 0.0;
      pTime[1] = pTime[0] + delta_time;


      // ??
      double SPD;
      int element_number = 0;
      int npt = 0; // integration point number

      // A. Create Properties vector
      const Properties & rMaterialProperties = rModelData.GetProperties();
      int number_properties;
      double* pPropertiesVector;
      this->CreateConstitutiveParametersVector( pPropertiesVector, number_properties, rMaterialProperties);


      // B. Create State Variables vector
      int number_state_variables =  this->GetNumberOfStateVariables();
      double* pStateVariables;
      this->CreateStateVariablesVector( pStateVariables, number_state_variables);

      // C. Create incremental strain vector
      double* pStrain;
      double* pDeltaStrain;
      this->CreateStrainsVectors( Variables, pStrain, pDeltaStrain);

      // D. Stress at the initial step
      double* pStressVector;
      this->CreateStressAtInitialState( Variables, pStressVector);


      // E. Set Umat constitutive law number
      int material_number = this->GetConstitutiveEquationNumber();


      bool DoingCorrection = false;
      if (this->mInitializedModel) {
         double PreviousMeanStress = 0;
         for (unsigned int i = 0; i < 3; i++)
            PreviousMeanStress += this->mStressVectorFinalized[i]/3.0;

         if (PreviousMeanStress > -2.5)
            DoingCorrection = true;
      }

      DoingCorrection = false;
      if ( DoingCorrection == false) {

         umat_wrapper_( pStressVector, pStateVariables, (double**) pConstitutiveMatrix, NULL, &SPD,
               NULL, NULL, NULL, NULL, NULL, pStrain, pDeltaStrain,
               pTime, &delta_time, NULL, NULL, NULL, NULL, NULL, &ndi, &nshr, &ntens, &number_state_variables, pPropertiesVector, &number_properties,
               NULL, NULL, NULL, NULL, NULL, NULL, &element_number, &npt, NULL, NULL, NULL, NULL, &material_number );
      }


      // Save stress vector
      VectorType StressVector = ZeroVector(6);
      Matrix Matrix = ZeroMatrix(6);

      if ( DoingCorrection) {

         double E(100);
         double nu(0.49);

         double diagonal =   E/(1.0+nu)/(1.0-2.0*nu) * (1.0-nu);
         double nodiagonal = E/(1.0+nu)/(1.0-2.0*nu) * ( nu);
         double corte      = E/(1.0+nu)/2.0;

         VectorType EpsVector(6);
         noalias( EpsVector) = ZeroVector(6);
         for (unsigned int i = 0; i < 6; i++)
            EpsVector(i) = pDeltaStrain[i];


         for (unsigned int i = 0; i<3; ++i) {
            for (unsigned int j = 0; j<3; ++j) {
               if (i == j) {
                  Matrix(i,i) = diagonal;
               }
               else {
                  Matrix(i,j) = nodiagonal;
               }
            }
         }

         for (unsigned int j = 3; j<6; ++j)
            Matrix(j,j) = corte;


         StressVector = this->mStressVectorFinalized + prod( Matrix , EpsVector);

         for (unsigned int i = 0; i < 6; i++)
            pStressVector[i] = StressVector(i);



      } else {

         for (unsigned int i = 0; i < 6; i++)
            StressVector(i) = pStressVector[i];


         // Save constitutive matrix
         for (unsigned int i = 0; i < 6; i++) {
            for (unsigned int j = 0; j< 6; j++) {
               Matrix(i,j) = pConstitutiveMatrix[i][j];
            }
         }

      }

      rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);
      this->SetConstitutiveMatrix( rConstitutiveMatrix, Matrix, rStressMatrix);


      // update internal variables
      if ( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES) ) {
         this->UpdateVariables(Variables, pStressVector, pStateVariables);
      }

      delete [] pPropertiesVector;
      delete [] pStateVariables;
      delete [] pStrain;
      delete [] pDeltaStrain;
      delete [] pStressVector;

      /*std::cout << " C " << rConstitutiveMatrix << std::endl;
        std::cout << " s " << rStressMatrix << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;*/

      KRATOS_CATCH(" ")
   }


   //********************************************************************************
   // Create the vector of constitutive parameters
   void SanisandSmallStrainUmatModel::CreateConstitutiveParametersVector(double* & pVector, int & rNumberParameters, const Properties & rMaterialProperties)
   {
      KRATOS_TRY

      rNumberParameters = 19;
      pVector = new double[rNumberParameters];

      pVector[0] = 100.0; // p_a atmospheric pressure
      pVector[1] = 0.8191; // e0
      pVector[2] = 0.00178; // lambda 
      pVector[3] = 2.4352; // epsi 
      pVector[4] = 1.287; // M_c or phi_c 
      pVector[5] = 1.0; // M_e or phi_e
      pVector[6] = 0.01; // m
      pVector[7] = 400.0; // Go
      pVector[8] = 0.05; // nu
      pVector[9] = 4.05; // h_0
      pVector[10] = 1.100; // c_h
      pVector[11] = 2.800; // n^b
      pVector[12] = 0.550; // A_0
      pVector[13] = 2.564; // n_d
      pVector[14] = 0; // z_max
      pVector[15] = 0; // c_z
      pVector[16] = 0; // K_w
      pVector[17] = 0.01; // p_tmult stabilising param
      pVector[18] = 0.6780; // e


      if ( rMaterialProperties.Has(P_ATM) )
      {
         if ( rMaterialProperties(P_ATM) > 0.0 )
            pVector[0] = rMaterialProperties[P_ATM];
      }
      if ( rMaterialProperties.Has(E0) )
      {
         if (rMaterialProperties[E0]> 0.0)
            pVector[1] = rMaterialProperties[E0];
      }
      if ( rMaterialProperties.Has(LAMBDA_SANISAND) )
      {
         if (rMaterialProperties[LAMBDA_SANISAND]> 0.0)
            pVector[2] = rMaterialProperties[LAMBDA_SANISAND];
      }
      if ( rMaterialProperties.Has(EPSI) )
      {
         if (rMaterialProperties[EPSI]> 0.0)
            pVector[3] = rMaterialProperties[EPSI];
      }
      if ( rMaterialProperties.Has(MC) )
      {
         if (rMaterialProperties[MC]> 0.0)
            pVector[4] = rMaterialProperties[MC];
      }
      if ( rMaterialProperties.Has(ME) )
      {
         if (rMaterialProperties[ME]> 0.0)
            pVector[5] = rMaterialProperties[ME];
      }
      if ( rMaterialProperties.Has(M_SANISAND) )
      {
         if (rMaterialProperties[M_SANISAND]> 0.0)
            pVector[6] = rMaterialProperties[M_SANISAND];
      }
      if ( rMaterialProperties.Has(G0_SANISAND) )
      {
         if (rMaterialProperties[G0_SANISAND]> 0.0)
            pVector[7] = rMaterialProperties[G0_SANISAND];
      }
      if ( rMaterialProperties.Has(NU_SANISAND) )
      {
         if (rMaterialProperties[NU_SANISAND]> 0.0)
            pVector[8] = rMaterialProperties[NU_SANISAND];
      }
      if ( rMaterialProperties.Has(H0) )
      {
         if (rMaterialProperties[H0]> 0.0)
            pVector[9] = rMaterialProperties[H0];
      }
      if ( rMaterialProperties.Has(CH) )
      {
         if (rMaterialProperties[CH]> 0.0)
            pVector[10] = rMaterialProperties[CH];
      }
      if ( rMaterialProperties.Has(NB) )
      {
         if (rMaterialProperties[NB]> 0.0)
            pVector[11] = rMaterialProperties[NB];
      }
      if ( rMaterialProperties.Has(A_SANISAND) )
      {
         if (rMaterialProperties[A_SANISAND]> 0.0)
            pVector[12] = rMaterialProperties[A_SANISAND];
      }
      if ( rMaterialProperties.Has(ND) )
      {
         if (rMaterialProperties[ND]> 0.0)
            pVector[13] = rMaterialProperties[ND];
      }
      if ( rMaterialProperties.Has(PTMULT) )
      {
         if (rMaterialProperties[PTMULT]> 0.0)
            pVector[17] = rMaterialProperties[PTMULT];
      }
      if ( rMaterialProperties.Has( VOID_RATIO) )
      {
         if (rMaterialProperties[VOID_RATIO]> 0.0)
            pVector[18] = rMaterialProperties[VOID_RATIO];
      }


      KRATOS_CATCH("")
   }

} // Namespace Kratos
