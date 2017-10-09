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
#include "custom_models/small_strain_Umat_model.hpp"

/* WRAPPER */
extern "C" void umat_wrapper_( double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
			       double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
			       double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
			       char* MATERL, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
			       double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
			       double** DFGRD1, double* NOEL, int* NPT, double* KSLAY, double* KSPT, double* KSTEP,
			       double* KINC, int* MATERIALNUMBER );



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SmallStrainUmatModel::SmallStrainUmatModel()
      : ConstitutiveModel()
   {
      mInitializedModel = false;
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SmallStrainUmatModel::SmallStrainUmatModel(const SmallStrainUmatModel& rOther)
      : ConstitutiveModel(rOther), mInitializedModel(rOther.mInitializedModel)
   {
         this->mInitializedModel = rOther.mInitializedModel; 
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer SmallStrainUmatModel::Clone() const
   {
      return ( SmallStrainUmatModel::Pointer(new SmallStrainUmatModel(*this)) );
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   SmallStrainUmatModel& SmallStrainUmatModel::operator=(SmallStrainUmatModel const& rOther)
   {
      ConstitutiveModel::operator=(rOther);
      this->mInitializedModel = rOther.mInitializedModel; 
      return *this;
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
         mpStateVariablesFinalized = new double[number_state_variables];
         for (int i = 0; i < number_state_variables; i++)
            mpStateVariablesFinalized[i] = 0.0;

         mpStressVectorFinalized = new double[6];
         mpStrainVectorFinalized = new double[6];
         for (unsigned int i = 0; i < 6; i++) {
            mpStressVectorFinalized[i] = 0.0;
            mpStrainVectorFinalized[i] = 0.0;
         }
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
   void SmallStrainUmatModel::InitializeElasticData(ModelDataType& rValues, ElasticDataType& rVariables)
   {
      KRATOS_TRY

      //set model data pointer
      rVariables.SetModelData(rValues);
      rVariables.SetState(rValues.State);     

      //add initial strain
      if(this->mOptions.Is(ConstitutiveModel::ADD_HISTORY_VECTOR) && this->mOptions.Is(ConstitutiveModel::HISTORY_STRAIN_MEASURE) ){
         VectorType StrainVector;
         StrainVector = ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);
         for(unsigned int i=0; i<StrainVector.size(); i++)
         {
            StrainVector[i] += this->mHistoryVector[i];	
         }
         rValues.StrainMatrix = ConstitutiveModelUtilities::StrainVectorToTensor(StrainVector, rValues.StrainMatrix);
      }

      rValues.SetStrainMeasure( ConstitutiveModelData::CauchyGreen_None);

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

      ElasticDataType Variables;
      this->InitializeElasticData(rValues,Variables);

      VectorType StrainVector;
      StrainVector = ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);

      // allocate variables
      int ndi = 3, nshr = 3, ntens = 6; // number of stress and strain components
      int npt = 0; // integration point number
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
      pTime[0] = 0.0;
      pTime[1] = pTime[0] + delta_time;
          

      // ??
      double SPD[1];

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
      this->CreateStrainsVectors( pStrain, pDeltaStrain, StrainVector);

      // D. Stress at the initial step
      double* pStressVector;
      this->CreateStressAtInitialState( pStressVector);


      // E. Set Umat constitutive law number
      int material_number = 0;


      umat_wrapper_( pStressVector, pStateVariables, (double**) pConstitutiveMatrix, NULL, SPD,
         NULL, NULL, NULL, NULL, NULL, pStrain, pDeltaStrain,
		   pTime, &delta_time, NULL, NULL, NULL, NULL, NULL, &ndi, &nshr, &ntens, &number_state_variables, pPropertiesVector, &number_properties,
		   NULL, NULL, NULL, NULL, NULL, NULL, NULL, &npt, NULL, NULL, NULL, NULL, &material_number );


      VectorType StressVector = ZeroVector(6);
      for (unsigned int i = 0; i < 6; i++)
         StressVector(i) = pStressVector[i];

      rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

      for (unsigned int i = 0; i < 6; i++)
         for (unsigned int j = 0; j < 6; j++)
            rConstitutiveMatrix(i,j) = pConstitutiveMatrix[i][j];

      if ( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES) )
         this->UpdateVariables( pStressVector, StrainVector, pStateVariables);


      KRATOS_CATCH(" ")
   }


   //************************************************************************************
   //************************************************************************************

   int SmallStrainUmatModel::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      return 0;

      KRATOS_CATCH(" ")
   }

} // Namespace Kratos
