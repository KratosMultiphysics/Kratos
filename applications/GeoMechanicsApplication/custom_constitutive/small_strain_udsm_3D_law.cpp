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

#include "custom_constitutive/small_strain_udsm_3D_law.hpp"
#include <algorithm>

#ifdef KRATOS_COMPILED_IN_WINDOWS
#include "windows.hpp"
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
#include <dlfcn.h>
#endif

namespace Kratos
{
/*
   - structure of the functions in PLAXIS UDSM:
   - Function to get stress, stiffness matrix, attribute, number of state variables, ...
   void User_Mod(int *IDTASK, int *IMOD, int *ISUNDR,
                 int *ISTEP, int *ITER, int *IEL, int *INT,
                 double *X, double *Y, double *Z,
                 double *TIME0, double *DTIME,
                 double *PROPS, double *SIG0, double *SWP0, double *STVAR0,
                 double *DEPS, double **D, double *BULKW,
                 double *SIG, double *SWP, double *STVAR, int *IPL,
                 int *NSTAT, int *NONSYM, int *ISTRSDEP, int *ITIMEDEP, int *ITANG,
                 int *IPRDIR, int *IPRJLEN, int *IABORT);

   void GetParamCount(int *IMOD, int *NPARAM);
   void GetStateVarCount(int *IMOD, int *NSTVAR);
*/

// calling convention (__cdecl, __stdcall, ...)
// __stdcall is the convention used by the WinAPI
#ifdef KRATOS_COMPILED_IN_WINDOWS
typedef void(__stdcall *f_GetParamCount)   (int *, int *);
typedef void(__stdcall *f_GetStateVarCount)(int *, int *);
typedef void(__stdcall *f_UserMod) (int    *, int     *, int    *,
                                    int    *, int     *, int    *, int *,
                                    double *, double  *, double *,
                                    double *, double  *,
                                    const double *, double  *, double *, double *,
                                    double *, double **, double *,
                                    double *, double  *, double *, int *,
                                    int    *, int     *, int    *, int *, int *,
                                    int    *, int     *, int    *);
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
typedef void(*f_GetParamCount)   (int *, int *);
typedef void(*f_GetStateVarCount)(int *, int *);
typedef void(*f_UserMod) (int    *, int     *, int    *,
                          int    *, int     *, int    *, int *,
                          double *, double  *, double *,
                          double *, double  *,
                          const double *, double  *, double *, double *,
                          double *, double **, double *,
                          double *, double  *, double *, int *,
                          int    *, int     *, int    *, int *, int *,
                          int    *, int     *, int    *);
#endif

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainUDSM3DLaw::SmallStrainUDSM3DLaw()
   : ConstitutiveLaw(),
     mIsModelInitialized(false),
     mIsUDSMLoaded(false)
   {
    KRATOS_TRY;

    KRATOS_CATCH("")

   }

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
SmallStrainUDSM3DLaw::SmallStrainUDSM3DLaw(const SmallStrainUDSM3DLaw &rOther)
   : ConstitutiveLaw(rOther),
     mStressVector(rOther.mStressVector),
     mStressVectorFinalized(rOther.mStressVectorFinalized),
     mDeltaStrainVector(rOther.mDeltaStrainVector),
     mStrainVectorFinalized(rOther.mStrainVectorFinalized),
     mIsModelInitialized(rOther.mIsModelInitialized),
     mIsUDSMLoaded(rOther.mIsUDSMLoaded),
     mAttributes(rOther.mAttributes),
     mStateVariables(rOther.mStateVariables),
     mStateVariablesFinalized(rOther.mStateVariablesFinalized)

{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::ConstitutiveLaw()") << std::endl;

   for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
      for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
         mMatrixD[i][j] = rOther.mMatrixD[i][j];

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::ConstitutiveLaw()") << std::endl;

   KRATOS_CATCH("");
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainUDSM3DLaw::Clone() const
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::Clone()") << std::endl;

   return Kratos::make_shared<SmallStrainUDSM3DLaw>(*this);

   KRATOS_CATCH("");
}

//********************************ASSIGNMENT******************************************
//************************************************************************************
SmallStrainUDSM3DLaw &SmallStrainUDSM3DLaw::operator=(SmallStrainUDSM3DLaw const &rOther)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::operator=()") << std::endl;

   ConstitutiveLaw::operator=(rOther);
   this->mIsModelInitialized      = rOther.mIsModelInitialized;
   this->mIsUDSMLoaded            = rOther.mIsUDSMLoaded;
   this->mAttributes              = rOther.mAttributes;
   this->mStateVariables          = rOther.mStateVariables;
   this->mStateVariablesFinalized = rOther.mStateVariablesFinalized;
   this->mStressVector            = rOther.mStressVector;
   this->mStressVectorFinalized   = rOther.mStressVectorFinalized;
   this->mDeltaStrainVector       = rOther.mDeltaStrainVector;
   this->mStrainVectorFinalized   = rOther.mStrainVectorFinalized;

   for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
      for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
         this->mMatrixD[i][j] = rOther.mMatrixD[i][j];

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::operator=()") << std::endl;

   return *this;

   KRATOS_CATCH("");
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainUDSM3DLaw::~SmallStrainUDSM3DLaw() {}

//***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
//************************************************************************************

void SmallStrainUDSM3DLaw::GetLawFeatures(Features &rFeatures)
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::GetLawFeatures()") << std::endl;

   //Set the type of law
   rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
   rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);

   if (mIsModelInitialized) {
      if (mAttributes[IS_NON_SYMMETRIC] == 1) {
         rFeatures.mOptions.Set(ANISOTROPIC);
      } else {
         rFeatures.mOptions.Set(ISOTROPIC);
      }
   } else {
      rFeatures.mOptions.Set(ISOTROPIC);
   }

   //Set strain measure required by the consitutive law
   rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
   //rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

   //Set the spacedimension
   rFeatures.mSpaceDimension = WorkingSpaceDimension();

   //Set the strain size
   rFeatures.mStrainSize = GetStrainSize();

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::GetLawFeatures()") << std::endl;

   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
int SmallStrainUDSM3DLaw::Check(const Properties &rMaterialProperties,
                                const GeometryType &rElementGeometry,
                                const ProcessInfo &rCurrentProcessInfo)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::Check()") << std::endl;

   // Verify Properties variables
   if (rMaterialProperties.Has(UDSM_NAME) == false || rMaterialProperties[UDSM_NAME] == "")
      KRATOS_THROW_ERROR(std::invalid_argument, 
                         "UDSM_NAME has Key zero, is not defined or has an invalid value for property",
                         rMaterialProperties.Id())

   if (rMaterialProperties.Has(UDSM_NUMBER) == false || rMaterialProperties[UDSM_NUMBER] <= 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "UDSM_NUMBER has Key zero, is not defined or has an invalid value for property",
                         rMaterialProperties.Id())
   if (rMaterialProperties.Has(IS_FORTRAN_UDSM) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "IS_FORTRAN_UDSM has Key zero, is not defined or has an invalid value for property",
                         rMaterialProperties.Id())

   // load UDSM model
   if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

   if (!mIsUDSMLoaded)
   {
      KRATOS_THROW_ERROR(std::runtime_error, "cannot load the specified UDSM ", rMaterialProperties[UDSM_NAME]);
   }

   const int nUmatParametersSize = rMaterialProperties[UMAT_PARAMETERS].size();
   const int nParametersUDSM = GetNumberOfMaterialParametersFromUDSM(rMaterialProperties);
   if ( nUmatParametersSize != nParametersUDSM)
   {
      KRATOS_THROW_ERROR(std::runtime_error, "Number of parameters is wrong."
                                             " The UDSM gives " + std::to_string(nParametersUDSM)
                                             + " while size of UMAT_PARAMETERS is " + std::to_string(nUmatParametersSize),
                                             rMaterialProperties[UDSM_NAME]);
   }


   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::Check()") << std::endl;
   return 0;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::InitializeMaterial(const Properties &rMaterialProperties,
                                              const GeometryType &rElementGeometry,
                                              const Vector &rShapeFunctionsValues)

{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::InitializeMaterial()") << std::endl;

   // loading the model
   mIsUDSMLoaded = loadUDSM(rMaterialProperties);

   ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::InitializeMaterial()") << std::endl;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::ResetStateVariables(const Properties& rMaterialProperties)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::ResetStateVariables()") << std::endl;

   // reset state variables
   int nStateVariables = GetNumberOfStateVariablesFromUDSM(rMaterialProperties);
   nStateVariables = std::max(nStateVariables, 1);
   mStateVariables.resize(nStateVariables);
   noalias(mStateVariables) = ZeroVector(nStateVariables);

   mStateVariablesFinalized.resize(nStateVariables);
   noalias(mStateVariablesFinalized) = ZeroVector(nStateVariables);

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::ResetStateVariables()") << std::endl;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                         const GeometryType& rElementGeometry,
                                         const Vector& rShapeFunctionsValues)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::ResetMaterial()") << std::endl;

   // reset state variables
   SetAttributes(rMaterialProperties);
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

   // state variables
   noalias(mStateVariables)         = ZeroVector(mStateVariables.size());
   noalias(mStateVariablesFinalized)= ZeroVector(mStateVariablesFinalized.size());

   mIsModelInitialized = false;

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::ResetMaterial()") << std::endl;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::SetAttributes(const Properties& rMaterialProperties)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::SetAttributes()") << std::endl;

   if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

   int IDTask = ATTRIBUTES;

   // process data
   double deltaTime = 0.0;
   double time      = 0.0;
   int    iStep     = 0;
   int    iteration = 0;

   // number of the model in the shared libaray (DLL)
   int modelNumber = rMaterialProperties[UDSM_NUMBER];

   // not needed:
   double bulkWater = 0.0;
   double excessPorePressurePrevious = 0.0;
   double excessPorePressureCurrent = 0.0;
   double X(0.0), Y(0.0), Z(0.0);
   int iElement = 0;
   int integrationNumber = 0;
   int iPlastic = 0;
   int isUndr = 0;
   int nStateVariables = 0;

   // variable to check if an error happend in the model:
   int iAbort = 0;
   int nSizeProjectDirectory = mProjectDirectory.size();
   std::vector<double> StateVariablesFinalized;
   std::vector<double> StateVariables;

   const auto &MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
   pUserMod(&IDTask, &modelNumber, &isUndr,
            &iStep, &iteration, &iElement, &integrationNumber,
            &X, &Y, &Z,
            &time, &deltaTime,
            &(MaterialParameters.data()[0]), &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious, 
            StateVariablesFinalized.data(),
            &(mDeltaStrainVector.data()[0]), (double **)mMatrixD, &bulkWater,
            &(mStressVector.data()[0]), &excessPorePressureCurrent, StateVariables.data(), &iPlastic,
            &nStateVariables, 
            &mAttributes[IS_NON_SYMMETRIC], &mAttributes[IS_STRESS_DEPENDENT],
            &mAttributes[IS_TIME_DEPENDENT], &mAttributes[USE_TANGENT_MATRIX],
            mProjectDirectory.data(), &nSizeProjectDirectory, 
            &iAbort);

   if (iAbort != 0)
   {
      // KRATOS_INFO("GetNumberOfStateVariablesFromUDSM, iAbort !=0")<< std::endl;
      KRATOS_THROW_ERROR(std::runtime_error, 
                         "the specified UDSM returns an error while call UDSM with IDTASK" + std::to_string(IDTask) + ". UDSM",
                         rMaterialProperties[UDSM_NAME]);
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::SetAttributes()") << std::endl;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
int SmallStrainUDSM3DLaw::GetNumberOfStateVariablesFromUDSM(const Properties& rMaterialProperties)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::GetNumberOfStateVariablesFromUDSM()") << std::endl;

   if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

   int IDTask = NUMBER_OF_STATE_VARIABLES;

   // process data
   double deltaTime = 0.0;
   double time      = 0.0;
   int    iStep     = 0;
   int    iteration = 0;

   // number of the model in the shared libaray (DLL)
   int modelNumber = rMaterialProperties[UDSM_NUMBER];

   // not needed:
   double bulkWater = 0.0;
   double excessPorePressurePrevious = 0.0;
   double excessPorePressureCurrent = 0.0;
   double X(0.0), Y(0.0), Z(0.0);
   int iElement = 0;
   int integrationNumber = 0;
   int iPlastic = 0;
   int isUndr = 0;
   int nStateVariables = 0;

   // variable to check if an error happend in the model:
   int iAbort = 0;
   int nSizeProjectDirectory = mProjectDirectory.size();
   std::vector<double> StateVariablesFinalized;
   std::vector<double> StateVariables;

   const auto &MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
   pUserMod(&IDTask, &modelNumber, &isUndr,
            &iStep, &iteration, &iElement, &integrationNumber,
            &X, &Y, &Z,
            &time, &deltaTime,
            &(MaterialParameters.data()[0]), &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious, 
            StateVariablesFinalized.data(),
            &(mDeltaStrainVector.data()[0]), (double **)mMatrixD, &bulkWater,
            &(mStressVector.data()[0]), &excessPorePressureCurrent, StateVariables.data(), &iPlastic,
            &nStateVariables, 
            &mAttributes[IS_NON_SYMMETRIC], &mAttributes[IS_STRESS_DEPENDENT],
            &mAttributes[IS_TIME_DEPENDENT], &mAttributes[USE_TANGENT_MATRIX],
            mProjectDirectory.data(), &nSizeProjectDirectory, 
            &iAbort);

   if (iAbort != 0)
   {
      // KRATOS_INFO("GetNumberOfStateVariablesFromUDSM, iAbort !=0")<< std::endl;
      KRATOS_THROW_ERROR(std::runtime_error, 
                         "the specified UDSM returns an error while call UDSM with IDTASK" + std::to_string(IDTask) + ". UDSM",
                         rMaterialProperties[UDSM_NAME]);
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::GetNumberOfStateVariablesFromUDSM()") << std::endl;

   return nStateVariables;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
int SmallStrainUDSM3DLaw::GetNumberOfMaterialParametersFromUDSM(const Properties& rMaterialProperties)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::GetNumberOfMaterialParametersFromUDSM()") << std::endl;

   if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

   int nUDSM = rMaterialProperties[UDSM_NUMBER];
   int nParameters(0);
   pGetParamCount(&nUDSM, &nParameters);

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::GetNumberOfMaterialParametersFromUDSM()") << std::endl;

   return nParameters;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
bool SmallStrainUDSM3DLaw::loadUDSM(const Properties &rMaterialProperties)
{
   KRATOS_TRY;
   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::loadUDSM()") << std::endl;

   bool isLoaded = false;

#ifdef KRATOS_COMPILED_IN_WINDOWS
   isLoaded = loadUDSMWindows(rMaterialProperties);
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::loadUDSM()") << std::endl;

   return isLoaded;
#endif

#if defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
   isLoaded = loadUDSMLinux(rMaterialProperties);
   return isLoaded;
#endif

   KRATOS_THROW_ERROR(std::logic_error, "loadUDSM is not supported yet for Mac OS applications", "");

   return isLoaded;

   KRATOS_CATCH(" ");
}

//----------------------------------------------------------------------------------------
bool SmallStrainUDSM3DLaw::loadUDSMLinux(const Properties &rMaterialProperties)
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
      KRATOS_INFO("Error in loadUDSMLinux") << "cannot load the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
      KRATOS_THROW_ERROR(std::runtime_error, "cannot load the specified UDSM ", rMaterialProperties[UDSM_NAME]);
      return false;
   }
   
   // resolve function GetParamCount address
   pGetParamCount = (f_GetParamCount)dlsym(lib_handle, "getparamcount");
   if (!pGetParamCount)
   {
      pGetParamCount = (f_GetParamCount)dlsym(lib_handle, "getparamcount_");
      if (!pGetParamCount)
      {
         KRATOS_INFO("Error in loadUDSMLinux") << "cannot load function GetParamCount in the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
         KRATOS_THROW_ERROR(std::runtime_error, "cannot load function GetParamCount in the specified UDSM ", rMaterialProperties[UDSM_NAME]);
         return false;
      }
   }

   // resolve function GetStateVarCount address
   pGetStateVarCount = (f_GetStateVarCount)dlsym(lib_handle, "getstatevarcount");

   pUserMod = (f_UserMod)dlsym(lib_handle, "user_mod");
   if (!pUserMod)
   {
      pUserMod = (f_UserMod)dlsym(lib_handle, "user_mod_");
      if (!pUserMod)
      {
         KRATOS_INFO("Error in loadUDSMLinux") << "cannot load function User_Mod in the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
         KRATOS_THROW_ERROR(std::runtime_error, "cannot load function User_Mod in the specified UDSM ", rMaterialProperties[UDSM_NAME]);
         return false;
      }
   }

   return true;

#else
   KRATOS_THROW_ERROR(std::logic_error, "loadUDSMLinux should be called in Linux applications", "");
   return false;
#endif
}

//----------------------------------------------------------------------------------------
bool SmallStrainUDSM3DLaw::loadUDSMWindows(const Properties &rMaterialProperties)
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::loadUDSMWindows()") << std::endl;

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
      KRATOS_INFO("Error in loadUDSMWindows") << "cannot load the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
      KRATOS_THROW_ERROR(std::runtime_error, "cannot load the specified UDSM ", rMaterialProperties[UDSM_NAME]);
      return false;
   }

   // resolve function GetParamCount address
   pGetParamCount = (f_GetParamCount)GetProcAddress(hGetProcIDDLL, "getparamcount");
   if (!pGetParamCount)
   {
      // check if the dll is compiled with gfortran
      pGetParamCount = (f_GetParamCount)GetProcAddress(hGetProcIDDLL, "getparamcount_");
      if (!pGetParamCount)
      {
         KRATOS_INFO("Error in loadUDSMWindows") << "cannot load function GetParamCount in the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
         KRATOS_THROW_ERROR(std::runtime_error, "cannot load function GetParamCount in the specified UDSM ", rMaterialProperties[UDSM_NAME]);
         return false;
      }
   }

   // resolve function GetStateVarCount address
   pGetStateVarCount = (f_GetStateVarCount)GetProcAddress(hGetProcIDDLL, "getstatevarcount");

   pUserMod = (f_UserMod)GetProcAddress(hGetProcIDDLL, "user_mod");
   if (!pUserMod)
   {
      // check if the dll is compiled with gfortran
      pUserMod = (f_UserMod)GetProcAddress(hGetProcIDDLL, "user_mod_");
      if (!pUserMod)
      {
         KRATOS_INFO("Error in loadUDSMWindows") << "cannot load function User_Mod in the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
         KRATOS_THROW_ERROR(std::runtime_error, "cannot load function User_Mod in the specified UDSM ", rMaterialProperties[UDSM_NAME]);
         return false;
      }
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::loadUDSMWindows()") << std::endl;

   return true;
#else
   KRATOS_THROW_ERROR(std::logic_error, "loadUDSMWindows should be called in Windows applications", "");
   return false;
#endif

   KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SmallStrainUDSM3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters & rValues)
{
   KRATOS_TRY;

   CalculateMaterialResponseCauchy(rValues);

   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters & rValues)
{
   KRATOS_TRY;

   CalculateMaterialResponseCauchy(rValues);

   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters & rValues)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CalculateMaterialResponseKirchhoff()") << std::endl;

   CalculateMaterialResponseCauchy(rValues);

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CalculateMaterialResponseKirchhoff()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CalculateMaterialResponseCauchy()") << std::endl;

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

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CalculateMaterialResponseCauchy()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues)
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::UpdateInternalDeltaStrainVector()") << std::endl;

   const Vector& rStrainVector = rValues.GetStrainVector();

   for (unsigned int i=0; i < mDeltaStrainVector.size(); ++i)
   {
      mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::UpdateInternalDeltaStrainVector()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::SetExternalStressVector(Vector& rStressVector)
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::SetExternalStressVector()") << std::endl;

   for (unsigned int i=0; i < rStressVector.size(); ++i)
   {
      rStressVector(i) = mStressVector[i];
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::SetExternalStressVector()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::SetInternalStressVector(const Vector& rStressVector)
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::SetInternalStressVector()") << std::endl;

   for (unsigned int i=0; i < mStressVectorFinalized.size(); ++i)
   {
      mStressVectorFinalized[i] = rStressVector(i);
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::SetInternalStressVector()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::SetInternalStrainVector()") << std::endl;

   for (unsigned int i=0; i < mStrainVectorFinalized.size(); ++i)
   {
      mStrainVectorFinalized[i] = rStrainVector(i);
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::SetInternalStrainVector()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CopyConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                   Matrix& rConstitutiveMatrix )
{
   KRATOS_TRY
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CopyConstitutiveMatrix()") << std::endl;

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

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CopyConstitutiveMatrix()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CalculateConstitutiveMatrix( ConstitutiveLaw::Parameters &rValues,
                                                        Matrix& rConstitutiveMatrix )
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CalculateConstitutiveMatrix()") << std::endl;

   // update strain vector
   UpdateInternalDeltaStrainVector(rValues);

   int IDTask = MATRIX_ELASTO_PLASTIC;

   CallUDSM(&IDTask, rValues);

   CopyConstitutiveMatrix(rValues, rConstitutiveMatrix);

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CalculateConstitutiveMatrix()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CalculateStress( ConstitutiveLaw::Parameters &rValues,
                                            Vector& rStressVector )
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CalculateStress()") << std::endl;

   // update strain vector
   UpdateInternalDeltaStrainVector(rValues);

   int IDTask = STRESS_CALCULATION;


   CallUDSM(&IDTask, rValues);

   SetExternalStressVector(rStressVector);

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CalculateStress()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CallUDSM(int *IDTask, ConstitutiveLaw::Parameters &rValues)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CallUDSM()") << std::endl;

   // process data
   double deltaTime = rValues.GetProcessInfo()[DELTA_TIME];
   double time      = rValues.GetProcessInfo()[TIME] - deltaTime;
   int    iStep     = rValues.GetProcessInfo()[STEP];
   int    iteration = rValues.GetProcessInfo()[NL_ITERATION_NUMBER];

   // number of the model in the shared libaray (DLL)
   const Properties& rMaterialProperties = rValues.GetMaterialProperties();
   int modelNumber = rMaterialProperties[UDSM_NUMBER];

   // number of state variables
   int nStateVariables = mStateVariablesFinalized.size();

   // not needed:
   double bulkWater = 0.0;
   double excessPorePressurePrevious = 0.0;
   double excessPorePressureCurrent = 0.0;
   double X(0.0), Y(0.0), Z(0.0);
   int iElement = 0;
   int integrationNumber = 0;
   int iPlastic = 0;
   int isUndr = 0;

   // variable to check if an error happend in the model:
   int iAbort = 0;
   int nSizeProjectDirectory = mProjectDirectory.size();

   const auto &MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
   pUserMod(IDTask, &modelNumber, &isUndr,
            &iStep, &iteration, &iElement, &integrationNumber,
            &X, &Y, &Z,
            &time, &deltaTime,
            &(MaterialParameters.data()[0]), &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious, 
            &(mStateVariablesFinalized.data()[0]),
            &(mDeltaStrainVector.data()[0]), (double **) mMatrixD, &bulkWater,
            &(mStressVector.data()[0]), &excessPorePressureCurrent, &(mStateVariables.data()[0]), &iPlastic,
            &nStateVariables, 
            &mAttributes[IS_NON_SYMMETRIC], &mAttributes[IS_STRESS_DEPENDENT],
            &mAttributes[IS_TIME_DEPENDENT], &mAttributes[USE_TANGENT_MATRIX],
            mProjectDirectory.data(), &nSizeProjectDirectory, 
            &iAbort);
   if (iAbort != 0)
   {
      KRATOS_INFO("CallUDSM, iAbort !=0")
                  << " iAbort: " << iAbort
                  << " the specified UDSM returns an error while call UDSM with IDTASK: " 
                  << std::to_string(*IDTask) << "." 
                  << " UDSM: " << rMaterialProperties[UDSM_NAME] 
                  << " UDSM_NUMBER: " << rMaterialProperties[UDSM_NUMBER]
                  << " Parameters: " << MaterialParameters
                  << std::endl;
      KRATOS_THROW_ERROR(std::runtime_error, 
                        "the specified UDSM returns an error while call UDSM with IDTASK: " 
                        + std::to_string(*IDTask) + ". UDSM: ",
                        rMaterialProperties[UDSM_NAME]);
   }


   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CallUDSM()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CallUDSM(int *IDTask, const Properties& rMaterialProperties)
{
   KRATOS_TRY;
   // KRATOS_INFO("01-SmallStrainUDSM3DLaw::CallUDSM()") << std::endl;

   // process data
   double deltaTime = 0.0;
   double time      = 0.0;
   int    iStep     = 0;
   int    iteration = 0;

   // number of the model in the shared libaray (DLL)
   int modelNumber = rMaterialProperties[UDSM_NUMBER];

   // number of state variables
   int nStateVariables = mStateVariablesFinalized.size();

   // not needed:
   double bulkWater = 0.0;
   double excessPorePressurePrevious = 0.0;
   double excessPorePressureCurrent = 0.0;
   double X(0.0), Y(0.0), Z(0.0);
   int iElement = 0;
   int integrationNumber = 0;
   int iPlastic = 0;
   int isUndr = 0;

   // variable to check if an error happend in the model:
   int iAbort = 0;
   int nSizeProjectDirectory = mProjectDirectory.size();
   
   const auto &MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
   pUserMod(IDTask, &modelNumber, &isUndr,
            &iStep, &iteration, &iElement, &integrationNumber,
            &X, &Y, &Z,
            &time, &deltaTime,
            &(MaterialParameters.data()[0]), &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious, 
            &(mStateVariablesFinalized.data()[0]),
            &(mDeltaStrainVector.data()[0]), (double **)mMatrixD, &bulkWater,
            &(mStressVector.data()[0]), &excessPorePressureCurrent, &(mStateVariables.data()[0]), &iPlastic,
            &nStateVariables,
            &mAttributes[IS_NON_SYMMETRIC], &mAttributes[IS_STRESS_DEPENDENT],
            &mAttributes[IS_TIME_DEPENDENT], &mAttributes[USE_TANGENT_MATRIX],
            mProjectDirectory.data(), &nSizeProjectDirectory, 
            &iAbort);

   if (iAbort != 0)
   {
      KRATOS_INFO("CallUDSM, iAbort !=0")
                  << " iAbort: " << iAbort
                  << " the specified UDSM returns an error while call UDSM with IDTASK: "
                  << std::to_string(*IDTask) << "." 
                  << " UDSM: " << rMaterialProperties[UDSM_NAME] 
                  << " UDSM_NUMBER: " << rMaterialProperties[UDSM_NUMBER]
                  << " Parameters: " << MaterialParameters
                  << std::endl;
      KRATOS_THROW_ERROR(std::runtime_error, 
                        "the specified UDSM returns an error while call UDSM with IDTASK" + std::to_string(*IDTask) + ". UDSM",
                        rMaterialProperties[UDSM_NAME]);
   }

   // KRATOS_INFO("11-SmallStrainUDSM3DLaw::CallUDSM()") << std::endl;
   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
   // Small deformation so we can call the Cauchy method
   InitializeMaterialResponseCauchy(rValues);

}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
   // Small deformation so we can call the Cauchy method
   InitializeMaterialResponseCauchy(rValues);
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
   // Small deformation so we can call the Cauchy method
   InitializeMaterialResponseCauchy(rValues);

}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
   KRATOS_TRY;
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::InitializeMaterialResponseCauchy()") << std::endl;

   if (!mIsModelInitialized)
   {
      // stress and strain vectors must be initialized:
      const Vector& rStressVector = rValues.GetStressVector();
      const Vector& rStrainVector = rValues.GetStrainVector();

      SetInternalStressVector(rStressVector);

      SetInternalStrainVector(rStrainVector);

      int IDTask = INITIALISATION;

      CallUDSM(&IDTask, rValues);

      mIsModelInitialized = true;
   }
   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::InitializeMaterialResponseCauchy()") << std::endl;

   KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters & rValues)
{
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::FinalizeMaterialResponseCauchy()") << std::endl;

   UpdateInternalStrainVectorFinalized(rValues);
   mStateVariablesFinalized = mStateVariables;
   mStressVectorFinalized   = mStressVector;

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::FinalizeMaterialResponseCauchy()") << std::endl;

}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters &rValues)
{
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::UpdateInternalStrainVectorFinalized()") << std::endl;
   const Vector& rStrainVector = rValues.GetStrainVector();

   for (unsigned int i=0; i < mStrainVectorFinalized.size(); ++i)
   {
      mStrainVectorFinalized[i] = rStrainVector(i);
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::UpdateInternalStrainVectorFinalized()") << std::endl;

}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::CalculateCauchyGreenStrain( ConstitutiveLaw::Parameters& rValues, 
                                                       Vector& rStrainVector )
{
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CalculateCauchyGreenStrain()") << std::endl;

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

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CalculateCauchyGreenStrain()") << std::endl;
}

//----------------------------------------------------------------------------------------
double& SmallStrainUDSM3DLaw::CalculateValue( ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<double>& rThisVariable,
                                              double& rValue )
{
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::CalculateValue()") << std::endl;

   Vector& rStrainVector = rParameterValues.GetStrainVector();
   Vector& rStressVector = rParameterValues.GetStressVector();

   if (rThisVariable == STRAIN_ENERGY)
   {
      this->CalculateCauchyGreenStrain(rParameterValues, rStrainVector);
      this->CalculateStress(rParameterValues, rStressVector);

      rValue = 0.5 * inner_prod( rStrainVector, rStressVector); // Strain energy = 0.5*E:C:E
   }

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::CalculateValue()") << std::endl;

   return( rValue );
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUDSM3DLaw::CalculateValue( ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<Vector>& rThisVariable,
                                              Vector& rValue )
{
   // KRATOS_INFO("01-SmallStrainUDSM3DLaw::CalculateValue()") << std::endl;

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
      CalculateMaterialResponseCauchy(rParameterValues);
      rValue = rParameterValues.GetStressVector();

      // Previous flags restored
      rFlags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flagConstTensor );
      rFlags.Set( ConstitutiveLaw::COMPUTE_STRESS, flagStress );
   }

   // KRATOS_INFO("11-SmallStrainUDSM3DLaw::CalculateValue()") << std::endl;

   return( rValue );
}

//----------------------------------------------------------------------------------------
Matrix& SmallStrainUDSM3DLaw::CalculateValue( ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<Matrix>& rThisVariable,
                                              Matrix& rValue )
{
   // KRATOS_INFO("02-SmallStrainUDSM3DLaw::CalculateValue()") << std::endl;

   if (rThisVariable == CONSTITUTIVE_MATRIX ||
       rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
       rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) 
   {
      this->CalculateConstitutiveMatrix(rParameterValues, rValue);
   }

   // KRATOS_INFO("12-SmallStrainUDSM3DLaw::CalculateValue()") << std::endl;

   return( rValue );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int SmallStrainUDSM3DLaw::GetStateVariableIndex(const Variable<double>& rThisVariable)
{
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::GetStateVariableIndex()") << std::endl;

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

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::GetStateVariableIndex()") << std::endl;

   return index - 1;
}

//----------------------------------------------------------------------------------------
Vector& SmallStrainUDSM3DLaw::GetValue( const Variable<Vector> &rThisVariable, Vector &rValue )
{
   // KRATOS_INFO("0-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

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

   // KRATOS_INFO("1-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

    return rValue;
}

//----------------------------------------------------------------------------------------
double& SmallStrainUDSM3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
   // KRATOS_INFO("01-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

   const int index = GetStateVariableIndex(rThisVariable);

   KRATOS_DEBUG_ERROR_IF( index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1) )
                        << "GetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

   rValue = mStateVariablesFinalized[index];

   // KRATOS_INFO("11-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

   return rValue;
}

//----------------------------------------------------------------------------------------
int& SmallStrainUDSM3DLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
   // KRATOS_INFO("02-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

   if (rThisVariable == NUMBER_OF_UMAT_STATE_VARIABLES)
   {
      rValue = mStateVariablesFinalized.size();
   }

   // KRATOS_INFO("12-SmallStrainUDSM3DLaw::GetValue()") << std::endl;

   return rValue;
}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::SetValue( const Variable<double>& rThisVariable,
                                     const double& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
   // KRATOS_INFO("01-SmallStrainUDSM3DLaw::SetValue()") << std::endl;

   const int index = GetStateVariableIndex(rThisVariable);

   KRATOS_DEBUG_ERROR_IF( index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1) )
                        << "SetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

   mStateVariablesFinalized[index] = rValue;

   // KRATOS_INFO("11-SmallStrainUDSM3DLaw::SetValue()") << std::endl;

}

//----------------------------------------------------------------------------------------
void SmallStrainUDSM3DLaw::SetValue( const Variable<Vector>& rThisVariable,
                                     const Vector& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
   // KRATOS_INFO("02-SmallStrainUDSM3DLaw::SetValue()") << std::endl;

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

   // KRATOS_INFO("12-SmallStrainUDSM3DLaw::SetValue()") << std::endl;

}

} // Namespace Kratos
