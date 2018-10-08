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
#include "custom_models/small_strain_umat_model.hpp"

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

   SmallStrainUmatModel::SmallStrainUmatModel()
      : ConstitutiveModel()
   {
      KRATOS_TRY

      mInitializedModel = false;

      KRATOS_CATCH("")
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SmallStrainUmatModel::SmallStrainUmatModel(const SmallStrainUmatModel& rOther)
      : ConstitutiveModel(rOther), mInitializedModel(rOther.mInitializedModel),
      mStateVariablesFinalized( rOther.mStateVariablesFinalized), mStressVectorFinalized( rOther.mStressVectorFinalized),
      mStrainVectorFinalized( rOther.mStrainVectorFinalized)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer SmallStrainUmatModel::Clone() const
   {
      KRATOS_TRY

          return Kratos::make_shared<SmallStrainUmatModel>(*this);

      KRATOS_CATCH("")
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   SmallStrainUmatModel& SmallStrainUmatModel::operator=(SmallStrainUmatModel const& rOther)
   {
      KRATOS_TRY

      ConstitutiveModel::operator=(rOther);
      this->mInitializedModel = rOther.mInitializedModel;
      return *this;

      KRATOS_CATCH("")
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SmallStrainUmatModel::~SmallStrainUmatModel()
   {
   }

   //***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
   //************************************************************************************

   void SmallStrainUmatModel::InitializeModel(ModelDataType& rValues)
   {
      KRATOS_TRY

      if ( mInitializedModel == false) {

         int number_state_variables = this->GetNumberOfStateVariables();
         mStateVariablesFinalized = ZeroVector(number_state_variables);

         mStressVectorFinalized.clear();
         mStrainVectorFinalized.clear();
      }
      mInitializedModel = true;

      KRATOS_CATCH(" ")
   }

   //************************************************************************************
   //************************************************************************************

   void SmallStrainUmatModel::FinalizeModel(ModelDataType& rValues)
   {
      KRATOS_TRY

      KRATOS_CATCH(" ")
   }

   //************************************************************************************
   //************************************************************************************
   void SmallStrainUmatModel::InitializeElasticData(ModelDataType& rValues, UmatDataType& rVariables)
   {
      KRATOS_TRY

      //set model data pointer
      rVariables.SetModelData(rValues);
      rVariables.SetState(rValues.State);

      //add initial strain
      if(this->mOptions.Is(ConstitutiveModel::ADD_HISTORY_VECTOR) && this->mOptions.Is(ConstitutiveModel::HISTORY_STRAIN_MEASURE) ){
         VectorType StrainVector;
         ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);
         for(unsigned int i=0; i<StrainVector.size(); i++)
         {
            StrainVector[i] += this->mHistoryVector[i];
         }
         rValues.StrainMatrix = ConstitutiveModelUtilities::StrainVectorToTensor(StrainVector, rValues.StrainMatrix);
      }

      rValues.SetStrainMeasure( ConstitutiveModelData::StrainMeasureType::CauchyGreen_None);

      rVariables.TotalStrainMatrix = rValues.StrainMatrix;
      rVariables.IncrementalDeformation = rValues.GetDeltaDeformationMatrix();

      KRATOS_CATCH(" ")
   }

   //************************************************************************************
   //************************************************************************************

   void SmallStrainUmatModel::CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
   {
      KRATOS_TRY

      KRATOS_ERROR << "calling the base class function in SmallStrainUmatModel ... illegal operation" << std::endl;

      KRATOS_CATCH(" ")
   }



   //************************************************************************************
   //************************************************************************************

   void SmallStrainUmatModel::CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
   {

      KRATOS_TRY
      Matrix ConstitutiveMatrix = ZeroMatrix(6);
      this->CalculateStressAndConstitutiveTensors( rValues, rStressMatrix, ConstitutiveMatrix);

      KRATOS_CATCH(" ")
   }


   //************************************************************************************
   //************************************************************************************

   void SmallStrainUmatModel::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
   {
      KRATOS_TRY

      MatrixType StressMatrix;
      this->CalculateStressAndConstitutiveTensors( rValues, StressMatrix, rConstitutiveMatrix);

      KRATOS_CATCH(" ")
   }


   //************************************************************************************
   //************************************************************************************

   void SmallStrainUmatModel::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
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
      const Properties & rMaterialProperties = rModelData.GetMaterialProperties();
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

      umat_wrapper_( pStressVector, pStateVariables, (double**) pConstitutiveMatrix, NULL, &SPD,
         NULL, NULL, NULL, NULL, NULL, pStrain, pDeltaStrain,
		   pTime, &delta_time, NULL, NULL, NULL, NULL, NULL, &ndi, &nshr, &ntens, &number_state_variables, pPropertiesVector, &number_properties,
		   NULL, NULL, NULL, NULL, NULL, NULL, &element_number, &npt, NULL, NULL, NULL, NULL, &material_number );


      // Save stress vector
      VectorType StressVector = ZeroVector(6);
      for (unsigned int i = 0; i < 6; i++)
         StressVector(i) = pStressVector[i];

      rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

      // Save constitutive matrix
      Matrix Matrix = ZeroMatrix(6);
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j< 6; j++) {
            Matrix(i,j) = pConstitutiveMatrix[i][j];
         }
      }

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

      KRATOS_CATCH(" ")
   }


   //************************************************************************************
   //************************************************************************************
   void SmallStrainUmatModel::SetConstitutiveMatrix( Matrix & rC, const Matrix & rpCBig, const MatrixType & rStressMatrix)

   {
      KRATOS_TRY

      if ( rC.size1() == 6 ) {
         for (unsigned int i = 0; i < 6; i++)
            for (unsigned int j = 0; j < 6; j++)
               rC(i,j) = rpCBig(i,j);

      } else if ( rC.size1() == 3) {
         rC(0,0) = rpCBig(0,0); rC(0,1) = rpCBig(0,1); rC(0,2) = rpCBig(0,3);
         rC(1,0) = rpCBig(1,0); rC(1,1) = rpCBig(1,1); rC(1,2) = rpCBig(1,3);
         rC(2,0) = rpCBig(3,0); rC(2,1) = rpCBig(3,1); rC(2,2) = rpCBig(3,3);

      } else {
         for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
               rC(i,j) = rpCBig(i,j);
            }
         }
      }


      KRATOS_CATCH("")
   }
   int SmallStrainUmatModel::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      return 0;

      KRATOS_CATCH(" ")
   }

} // Namespace Kratos
