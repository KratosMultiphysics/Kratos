// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes

// External includes

#include "custom_constitutive/small_strain_umat_3D_law.hpp"

#ifdef KRATOS_COMPILED_IN_WINDOWS
#include "windows.hpp"
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
#include <dlfcn.h>
#endif


namespace Kratos
{

#ifdef KRATOS_COMPILED_IN_WINDOWS
typedef void(__stdcall *f_UMATMod) (double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                                    double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                                    double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                                    char* MATERL, int* NDI, int* NSHR, int* NTENS, int* NSTATV, const double* PROPS, int* NPROPS,
                                    double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
                                    double** DFGRD1, int* NOEL, int* NPT, double* KSLAY, double* KSPT, int* KSTEP,
                                    int* KINC);
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
typedef void(*f_UMATMod) (double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                          double* rpl, double* ddsddt, double* drplde, double* drpldt, double* stran, double* dstran,
                          double* time, double* dtime, double* temp, double* dtemp, double* predef, double* dpred,
                          char* materl, int* ndi, int* nshr, int* ntens, int* nstatv, const double* props, int* nprops,
                          double* coords, double** drot, double* pnewdt, double* celent, double** dfgrd0,
                          double** dfgrd1, int* noel, int* npt, double* kslay, double* kspt, int* kstep,
                          int* kinc);
#endif

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT3DLaw::SmallStrainUMAT3DLaw()
   : ConstitutiveLaw(),
     mIsModelInitialized(false),
     mIsUMATLoaded(false)
   {
    KRATOS_TRY;

    KRATOS_CATCH("")

   }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
SmallStrainUMAT3DLaw::SmallStrainUMAT3DLaw(const SmallStrainUMAT3DLaw &rOther)
   : ConstitutiveLaw(rOther),
     mStressVector(rOther.mStressVector),
     mStressVectorFinalized(rOther.mStressVectorFinalized),
     mDeltaStrainVector(rOther.mDeltaStrainVector),
     mStrainVectorFinalized(rOther.mStrainVectorFinalized),
     mIsModelInitialized(rOther.mIsModelInitialized),
     mIsUMATLoaded(rOther.mIsUMATLoaded),
     mStateVariables(rOther.mStateVariables),
     mStateVariablesFinalized(rOther.mStateVariablesFinalized)

{
   KRATOS_TRY;

   for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
      for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
         mMatrixD[i][j] = rOther.mMatrixD[i][j];

   KRATOS_CATCH("");
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainUMAT3DLaw::Clone() const
{
   KRATOS_TRY;

   return Kratos::make_shared<SmallStrainUMAT3DLaw>(*this);

   KRATOS_CATCH("");
}

//********************************ASSIGNMENT******************************************
//************************************************************************************
SmallStrainUMAT3DLaw &SmallStrainUMAT3DLaw::operator=(SmallStrainUMAT3DLaw const &rOther)
{
   KRATOS_TRY;

   ConstitutiveLaw::operator=(rOther);
   this->mIsModelInitialized      = rOther.mIsModelInitialized;
   this->mIsUMATLoaded            = rOther.mIsUMATLoaded;
   this->mStateVariables          = rOther.mStateVariables;
   this->mStateVariablesFinalized = rOther.mStateVariablesFinalized;
   this->mStressVector            = rOther.mStressVector;
   this->mStressVectorFinalized   = rOther.mStressVectorFinalized;
   this->mDeltaStrainVector       = rOther.mDeltaStrainVector;
   this->mStrainVectorFinalized   = rOther.mStrainVectorFinalized;

   for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
      for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
         this->mMatrixD[i][j] = rOther.mMatrixD[i][j];

   return *this;

   KRATOS_CATCH("");
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainUMAT3DLaw::~SmallStrainUMAT3DLaw() {}

//***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
//************************************************************************************

void SmallStrainUMAT3DLaw::GetLawFeatures(Features &rFeatures)
{
   //Set the type of law
   rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
   rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);

   rFeatures.mOptions.Set(ISOTROPIC);


   //Set strain measure required by the consitutive law
   rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
   //rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

   //Set the spacedimension
   rFeatures.mSpaceDimension = WorkingSpaceDimension();

   //Set the strain size
   rFeatures.mStrainSize = GetStrainSize();
}

int SmallStrainUMAT3DLaw::Check(const Properties &rMaterialProperties,
                                const GeometryType &rElementGeometry,
                                const ProcessInfo &rCurrentProcessInfo)
{
   // Verify Properties variables
   if (rMaterialProperties.Has(UDSM_NAME) == false || rMaterialProperties[UDSM_NAME] == "")
      KRATOS_ERROR << "UDSM_NAME has Key zero, is not defined for property"
                   << rMaterialProperties.Id()
                   << std::endl;



   if (rMaterialProperties.Has(IS_FORTRAN_UDSM) == false)
      KRATOS_ERROR << "IS_FORTRAN_UDSM is not defined for property"
                   << rMaterialProperties.Id()
                   << std::endl;

   // load UMAT model
   if (!mIsUMATLoaded) mIsUMATLoaded = loadUMAT(rMaterialProperties);

   if (!mIsUMATLoaded)
   {
      KRATOS_ERROR << "cannot load the specified UMAT" << rMaterialProperties[UDSM_NAME] << std::endl;
   }

   return 0;
}

void SmallStrainUMAT3DLaw::InitializeMaterial(const Properties &rMaterialProperties,
                                              const GeometryType &rElementGeometry,
                                              const Vector &rShapeFunctionsValues)

{
   KRATOS_TRY;
   // we need to check if the model is loaded or not
   mIsUMATLoaded = loadUMAT(rMaterialProperties);

   ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

   KRATOS_CATCH(" ");
}

void SmallStrainUMAT3DLaw::ResetStateVariables(const Properties& rMaterialProperties)
{
   KRATOS_TRY;
   // reset state variables

   const auto &StateVariables = rMaterialProperties[STATE_VARIABLES];
   const unsigned int nStateVariables = StateVariables.size();

   mStateVariables.resize(nStateVariables);
   mStateVariablesFinalized.resize(nStateVariables);

   noalias(mStateVariables)          = StateVariables;
   noalias(mStateVariablesFinalized) = StateVariables;

   KRATOS_CATCH(" ");
}


void SmallStrainUMAT3DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                         const GeometryType& rElementGeometry,
                                         const Vector& rShapeFunctionsValues)
{
   KRATOS_TRY;

   // reset state variables
   ResetStateVariables(rMaterialProperties);

   // set stress vectors:
   noalias(mStressVector)          = ZeroVector(mStressVector.size());
   noalias(mStressVectorFinalized) = ZeroVector(mStressVectorFinalized.size());

   // set strain vectors:
   noalias(mDeltaStrainVector)     = ZeroVector(mDeltaStrainVector.size());
   noalias(mStrainVectorFinalized) = ZeroVector(mStrainVectorFinalized.size());

   for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
      for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
         mMatrixD[i][j] = 0.0;

   mIsModelInitialized = false;

   KRATOS_CATCH(" ");
}

bool SmallStrainUMAT3DLaw::loadUMAT(const Properties &rMaterialProperties)
{
   KRATOS_TRY;

   bool isLoaded = false;

#ifdef KRATOS_COMPILED_IN_WINDOWS
   return loadUMATWindows(rMaterialProperties);
#endif

#if defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
   return loadUMATLinux(rMaterialProperties);
#endif

   KRATOS_ERROR << "loadUMAT is not supported yet for Mac OS applications" << std::endl;

   return isLoaded;

   KRATOS_CATCH(" ");
}

bool SmallStrainUMAT3DLaw::loadUMATLinux(const Properties &rMaterialProperties)
{
#ifdef KRATOS_COMPILED_IN_LINUX
   void *lib_handle;

   lib_handle = dlopen((rMaterialProperties[UDSM_NAME]).c_str(), RTLD_LAZY);
   if (!lib_handle)
   {
      std::string name = rMaterialProperties[UDSM_NAME];
      // check if the name of the file is based on Windows extension
      std::size_t found = name.find(".dll");
      if (found!=std::string::npos)
      {
         // check if there is an equivalent .so file
         name.replace(found, 4, ".so");
         lib_handle = dlopen(name.c_str(), RTLD_LAZY);
      }
   }

   if (!lib_handle)
   {
      KRATOS_ERROR << "cannot load the specified UMAT " << rMaterialProperties[UDSM_NAME] << std::endl;
      return false;
   }
   if (rMaterialProperties[IS_FORTRAN_UDSM])
   {
      pUserMod = (f_UMATMod)dlsym(lib_handle, "umat_");
   }
   else
   {
      pUserMod = (f_UMATMod)dlsym(lib_handle, "umat");
   }


   if (!pUserMod)
   {
      KRATOS_ERROR << "cannot load function User_Mod in the specified UMAT " << rMaterialProperties[UDSM_NAME] << std::endl;
      return false;
   }

   return true;

#else
   KRATOS_ERROR << "loadUMATLinux should be called in Linux applications" << rMaterialProperties[UDSM_NAME] << std::endl;
   return false;
#endif
}

bool SmallStrainUMAT3DLaw::loadUMATWindows(const Properties &rMaterialProperties)
{
#ifdef KRATOS_COMPILED_IN_WINDOWS

   HINSTANCE hGetProcIDDLL = LoadLibrary((rMaterialProperties[UDSM_NAME]).c_str());

   if (!hGetProcIDDLL)
   {
      std::string name = rMaterialProperties[UDSM_NAME];
      // check if the name of the file is based on Linux extension
      std::size_t found = name.find(".so");
      if (found!=std::string::npos)
      {
         // check if there is an equivalent .dll file
         name.replace(found, 3, ".dll");
         hGetProcIDDLL = LoadLibrary(name.c_str());
      }
   }

   if (!hGetProcIDDLL)
   {
      KRATOS_INFO("Error in loadUMATWindows") << "cannot load the specified UMAT: " << rMaterialProperties[UDSM_NAME] << std::endl;
      KRATOS_ERROR << "cannot load the specified UMAT " << rMaterialProperties[UDSM_NAME] << std::endl;
      return false;
   }

   pUserMod = (f_UMATMod)GetProcAddress(hGetProcIDDLL, "umat");
   if (!pUserMod)
   {
      KRATOS_INFO("Error in loadUMATWindows") << "cannot load function umat in the specified UMAT: "<< rMaterialProperties[UDSM_NAME] << std::endl;
      KRATOS_ERROR << "cannot load function umat in the specified UMAT " << rMaterialProperties[UDSM_NAME] << std::endl;
      return false;
   }

   return true;
#else
   KRATOS_ERROR << "loadUMATWindows should be called in Windows applications" << rMaterialProperties[UDSM_NAME] << std::endl;
   return false;
#endif
}


//************************************************************************************
//************************************************************************************

void SmallStrainUMAT3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters & rValues)
{
   KRATOS_TRY;

   CalculateMaterialResponseCauchy(rValues);

   KRATOS_CATCH("");
}

void SmallStrainUMAT3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters & rValues)
{
   KRATOS_TRY;

   CalculateMaterialResponseCauchy(rValues);

   KRATOS_CATCH("");
}

void SmallStrainUMAT3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters & rValues)
{
   KRATOS_TRY;

   CalculateMaterialResponseCauchy(rValues);

   KRATOS_CATCH("");
}

void SmallStrainUMAT3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
   KRATOS_TRY;

   // Get Values to compute the constitutive law:
   Flags &rOptions=rValues.GetOptions();

   //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
   if (rOptions.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
   {
      Vector& rStrainVector = rValues.GetStrainVector();
      CalculateCauchyGreenStrain( rValues, rStrainVector);
   }

   if (rOptions.Is( ConstitutiveLaw::COMPUTE_STRESS ))
   {
      Vector& rStressVector = rValues.GetStressVector();
      CalculateStress(rValues, rStressVector);
   }

   if (rOptions.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) 
   {
      // Constitutive matrix (D matrix)
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
      CalculateConstitutiveMatrix(rValues, rConstitutiveMatrix);
   }

   KRATOS_CATCH("");
}


void SmallStrainUMAT3DLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();

   for (unsigned int i=0; i < mDeltaStrainVector.size(); ++i)
   {
      mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
   }
}

void SmallStrainUMAT3DLaw::SetExternalStressVector(Vector& rStressVector)
{
   for (unsigned int i=0; i < rStressVector.size(); ++i)
   {
      rStressVector(i) = mStressVector[i];
   }
}

void SmallStrainUMAT3DLaw::SetInternalStressVector(const Vector& rStressVector)
{
   for (unsigned int i=0; i < mStressVectorFinalized.size(); ++i)
   {
      mStressVectorFinalized[i] = rStressVector(i);
   }
}

void SmallStrainUMAT3DLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
   for (unsigned int i=0; i < mStrainVectorFinalized.size(); ++i)
   {
      mStrainVectorFinalized[i] = rStrainVector(i);
   }
}

void SmallStrainUMAT3DLaw::CopyConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                   Matrix& rConstitutiveMatrix )
{
   if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM])
   {
      // transfer fortran style matrix to C++ style
      for (unsigned int i = 0; i < VOIGT_SIZE_3D; i++) {
         for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[j][i];
         }
      }
   }
   else
   {
      for (unsigned int i = 0; i < VOIGT_SIZE_3D; i++) {
         for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) {
            rConstitutiveMatrix(i,j) = mMatrixD[i][j];
         }
      }
   }
}


void SmallStrainUMAT3DLaw::CalculateConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                        Matrix& rConstitutiveMatrix )
{
   KRATOS_TRY;

   // update strain vector
   UpdateInternalDeltaStrainVector(rValues);

   CallUMAT( rValues);

   CopyConstitutiveMatrix(rValues, rConstitutiveMatrix);

   KRATOS_CATCH("");
}

void SmallStrainUMAT3DLaw::CalculateStress( ConstitutiveLaw::Parameters &rValues,
                                            Vector& rStressVector )
{
   KRATOS_TRY;

   // update strain vector
   UpdateInternalDeltaStrainVector(rValues);

   CallUMAT(rValues);

   SetExternalStressVector(rStressVector);

   KRATOS_CATCH("");
}

void SmallStrainUMAT3DLaw::CallUMAT( ConstitutiveLaw::Parameters &rValues)
{
   KRATOS_TRY;

   // process data
   double deltaTime = rValues.GetProcessInfo()[DELTA_TIME];
   double time      = rValues.GetProcessInfo()[TIME] - deltaTime;
   int    iStep     = rValues.GetProcessInfo()[STEP];
   int    iteration = rValues.GetProcessInfo()[NL_ITERATION_NUMBER];

   // number of the model in the shared libaray (DLL)

   // number of state variables
   int nStateVariables = mStateVariablesFinalized.size();

   // not needed:
   int iElement = 0;
   int integrationNumber = 0;
   double SPD;// specific plastic dissipation
   double SSE;// ?
   double SCD;// ?
   char materialName;

   int ndi = N_DIM_3D;
   int nshr = 3;
   int ntens = VOIGT_SIZE_3D; // ??@@


   // stresses and state variables in the beginning of the steps needs to be given:
   mStressVector  = mStressVectorFinalized;
   mStateVariables = mStateVariablesFinalized;

   // variable to check if an error happend in the model:
   const auto &MaterialParameters = rValues.GetMaterialProperties()[UMAT_PARAMETERS];
   int nProperties = MaterialParameters.size();
   pUserMod(&(mStressVector.data()[0]), &(mStateVariables.data()[0]), (double **)mMatrixD,  &SSE,   &SPD,                          &SCD,
            NULL,                 NULL,                   NULL,                 NULL,   &(mStrainVectorFinalized.data()[0]), &(mDeltaStrainVector.data()[0]),
            &time,                &deltaTime,             NULL,                 NULL,   NULL,                          NULL,
            &materialName,        &ndi,                   &nshr,                &ntens, &nStateVariables,              &(MaterialParameters.data()[0]),
            &nProperties,         NULL,                   NULL,                 NULL,   NULL,                          NULL,
            NULL,                 &iElement,              &integrationNumber,   NULL,   NULL,                          &iStep,
            &iteration);

   KRATOS_CATCH("");

}

void SmallStrainUMAT3DLaw::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
   // Small deformation so we can call the Cauchy method
   InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUMAT3DLaw::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
   // Small deformation so we can call the Cauchy method
   InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUMAT3DLaw::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
   // Small deformation so we can call the Cauchy method
   InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUMAT3DLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
   KRATOS_TRY;

   if (!mIsModelInitialized)
   {
      // stress and strain vectors must be initialized:
      const Vector& rStressVector = rValues.GetStressVector();
      const Vector& rStrainVector = rValues.GetStrainVector();

      SetInternalStressVector(rStressVector);

      SetInternalStrainVector(rStrainVector);

      CallUMAT(rValues);
      mIsModelInitialized = true;
   }

   KRATOS_CATCH("");
}

void SmallStrainUMAT3DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUMAT3DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUMAT3DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUMAT3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters & rValues)
{
   UpdateInternalStrainVectorFinalized(rValues);
   mStateVariablesFinalized = mStateVariables;
   mStressVectorFinalized   = mStressVector;
}

void SmallStrainUMAT3DLaw::
   UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters &rValues)
{
   const Vector& rStrainVector = rValues.GetStrainVector();
   this->SetInternalStrainVector(rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainUMAT3DLaw::CalculateCauchyGreenStrain( ConstitutiveLaw::Parameters& rValues, 
                                                       Vector& rStrainVector )
{
   const SizeType space_dimension = this->WorkingSpaceDimension();

   //-Compute total deformation gradient
   const Matrix& F = rValues.GetDeformationGradientF();
   KRATOS_DEBUG_ERROR_IF(F.size1()!= space_dimension || F.size2() != space_dimension)
                         << "expected size of F " << space_dimension 
                         << "x" << space_dimension 
                         << ", got " << F.size1() 
                         << "x" << F.size2() << std::endl;

   Matrix E_tensor = prod(trans(F), F);
   for (unsigned int i=0; i<space_dimension; ++i)
      E_tensor(i,i) -= 1.0;
   E_tensor *= 0.5;

   noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}


double& SmallStrainUMAT3DLaw::CalculateValue( ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<double>& rThisVariable,
                                              double& rValue )
{
   Vector& rStrainVector = rParameterValues.GetStrainVector();
   Vector& rStressVector = rParameterValues.GetStressVector();

   if (rThisVariable == STRAIN_ENERGY) 
   {
      this->CalculateCauchyGreenStrain(rParameterValues, rStrainVector);
      this->CalculateStress(rParameterValues, rStressVector);

      rValue = 0.5 * inner_prod( rStrainVector, rStressVector); // Strain energy = 0.5*E:C:E
   }

   return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SmallStrainUMAT3DLaw::CalculateValue( ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<Vector>& rThisVariable,
                                              Vector& rValue )
{
   if (rThisVariable == STRAIN ||
       rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
       rThisVariable == ALMANSI_STRAIN_VECTOR) 
   {
      this->CalculateCauchyGreenStrain( rParameterValues, rValue);

   } else if (rThisVariable == STRESSES ||
              rThisVariable == CAUCHY_STRESS_VECTOR ||
              rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
              rThisVariable == PK2_STRESS_VECTOR) 
   {
        // Get Values to compute the constitutive law:
      Flags& rFlags = rParameterValues.GetOptions();

      // Previous flags saved
      const bool flagConstTensor = rFlags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
      const bool flagStress = rFlags.Is( ConstitutiveLaw::COMPUTE_STRESS );

      rFlags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
      rFlags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

      // We compute the stress
      SmallStrainUMAT3DLaw::CalculateMaterialResponseCauchy(rParameterValues);
      rValue = rParameterValues.GetStressVector();

      // Previous flags restored
      rFlags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flagConstTensor );
      rFlags.Set( ConstitutiveLaw::COMPUTE_STRESS, flagStress );
   }

   return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& SmallStrainUMAT3DLaw::CalculateValue( ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<Matrix>& rThisVariable,
                                              Matrix& rValue )
{
   if (rThisVariable == CONSTITUTIVE_MATRIX ||
       rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
       rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) 
   {
      this->CalculateConstitutiveMatrix(rParameterValues, rValue);
   }

   return( rValue );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int SmallStrainUMAT3DLaw::GetStateVariableIndex(const Variable<double>& rThisVariable)
{

   int index = -1;

   if (rThisVariable == STATE_VARIABLE_1)
       index = 1;
    else if (rThisVariable == STATE_VARIABLE_2)
       index = 2;
    else if (rThisVariable == STATE_VARIABLE_3)
       index = 3;
    else if (rThisVariable == STATE_VARIABLE_4)
       index = 4;
    else if (rThisVariable == STATE_VARIABLE_5)
       index = 5;
    else if (rThisVariable == STATE_VARIABLE_6)
       index = 6;
    else if (rThisVariable == STATE_VARIABLE_7)
       index = 7;
    else if (rThisVariable == STATE_VARIABLE_8)
       index = 8;
    else if (rThisVariable == STATE_VARIABLE_9)
       index = 9;

    else if (rThisVariable == STATE_VARIABLE_10)
       index = 10;
    else if (rThisVariable == STATE_VARIABLE_11)
       index = 11;
    else if (rThisVariable == STATE_VARIABLE_12)
       index = 12;
    else if (rThisVariable == STATE_VARIABLE_13)
       index = 13;
    else if (rThisVariable == STATE_VARIABLE_14)
       index = 14;
    else if (rThisVariable == STATE_VARIABLE_15)
       index = 15;
    else if (rThisVariable == STATE_VARIABLE_16)
       index = 16;
    else if (rThisVariable == STATE_VARIABLE_17)
       index = 17;
    else if (rThisVariable == STATE_VARIABLE_18)
       index = 18;
    else if (rThisVariable == STATE_VARIABLE_19)
       index = 19;
    else if (rThisVariable == STATE_VARIABLE_20)
       index = 20;

    else if (rThisVariable == STATE_VARIABLE_21)
       index = 21;
    else if (rThisVariable == STATE_VARIABLE_22)
       index = 22;
    else if (rThisVariable == STATE_VARIABLE_23)
       index = 23;
    else if (rThisVariable == STATE_VARIABLE_24)
       index = 24;
    else if (rThisVariable == STATE_VARIABLE_25)
       index = 25;
    else if (rThisVariable == STATE_VARIABLE_26)
       index = 26;
    else if (rThisVariable == STATE_VARIABLE_27)
       index = 27;
    else if (rThisVariable == STATE_VARIABLE_28)
       index = 28;
    else if (rThisVariable == STATE_VARIABLE_29)
       index = 29;

    else if (rThisVariable == STATE_VARIABLE_30)
       index = 30;
    else if (rThisVariable == STATE_VARIABLE_31)
       index = 31;
    else if (rThisVariable == STATE_VARIABLE_32)
       index = 32;
    else if (rThisVariable == STATE_VARIABLE_33)
       index = 33;
    else if (rThisVariable == STATE_VARIABLE_34)
       index = 34;
    else if (rThisVariable == STATE_VARIABLE_35)
       index = 35;
    else if (rThisVariable == STATE_VARIABLE_36)
       index = 36;
    else if (rThisVariable == STATE_VARIABLE_37)
       index = 37;
    else if (rThisVariable == STATE_VARIABLE_38)
       index = 38;
    else if (rThisVariable == STATE_VARIABLE_39)
       index = 39;

    else if (rThisVariable == STATE_VARIABLE_40)
       index = 40;
    else if (rThisVariable == STATE_VARIABLE_41)
       index = 41;
    else if (rThisVariable == STATE_VARIABLE_42)
       index = 42;
    else if (rThisVariable == STATE_VARIABLE_43)
       index = 43;
    else if (rThisVariable == STATE_VARIABLE_44)
       index = 44;
    else if (rThisVariable == STATE_VARIABLE_45)
       index = 45;
    else if (rThisVariable == STATE_VARIABLE_46)
       index = 46;
    else if (rThisVariable == STATE_VARIABLE_47)
       index = 47;
    else if (rThisVariable == STATE_VARIABLE_48)
       index = 48;
    else if (rThisVariable == STATE_VARIABLE_49)
       index = 49;

    else if (rThisVariable == STATE_VARIABLE_50)
       index = 50;

   return index -1;
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUMAT3DLaw::GetValue( const Variable<Vector> &rThisVariable, Vector &rValue )
{

   if (rThisVariable == STATE_VARIABLES)
   {
      if (rValue.size() != mStateVariablesFinalized.size())
         rValue.resize(mStateVariablesFinalized.size());

      noalias(rValue) = mStateVariablesFinalized;
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() != mStressVectorFinalized.size())
         rValue.resize(mStressVectorFinalized.size());

      noalias(rValue) = mStressVectorFinalized;
   }

    return rValue;
}

//----------------------------------------------------------------------------------------
double& SmallStrainUMAT3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{

   int index = GetStateVariableIndex(rThisVariable);

   KRATOS_DEBUG_ERROR_IF( index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1) )
                        << "GetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

   rValue = mStateVariablesFinalized[index];

    return rValue;
}

//----------------------------------------------------------------------------------------
int& SmallStrainUMAT3DLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if (rThisVariable == NUMBER_OF_UMAT_STATE_VARIABLES)
    {
       rValue = mStateVariablesFinalized.size();
    }

    return rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUMAT3DLaw::SetValue( const Variable<double>& rThisVariable,
                                     const double& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
   const int index = GetStateVariableIndex(rThisVariable);

   KRATOS_DEBUG_ERROR_IF( index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1) )
                        << "SetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

   mStateVariablesFinalized[index] = rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUMAT3DLaw::SetValue( const Variable<Vector>& rThisVariable,
                                     const Vector& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
   if (rThisVariable == STATE_VARIABLES)
   {
      if (rValue.size() == mStateVariablesFinalized.size()) 
      {
         for (unsigned int i=0; i < rValue.size(); ++i)
         {
            mStateVariablesFinalized[i] = rValue[i];
         }
      }
   }
   else if (rThisVariable == CAUCHY_STRESS_VECTOR)
   {
      if (rValue.size() == mStressVectorFinalized.size()) 
      {
         for (unsigned int i=0; i < rValue.size(); ++i)
         {
            mStressVectorFinalized[i] = rValue[i];
         }
      }
   }

}
} // Namespace Kratos
