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
#include <algorithm>

#include "custom_constitutive/small_strain_udsm_3D_law.hpp"
#include "custom_utilities/constitutive_law_utilities.hpp"

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
typedef void(__stdcall* f_GetParamCount)(int*, int*);
typedef void(__stdcall* f_GetStateVarCount)(int*, int*);
typedef void(__stdcall* f_UserMod)(int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   double*,
                                   double*,
                                   double*,
                                   double*,
                                   double*,
                                   const double*,
                                   double*,
                                   double*,
                                   double*,
                                   double*,
                                   double**,
                                   double*,
                                   double*,
                                   double*,
                                   double*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*,
                                   int*);
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
typedef void (*f_GetParamCount)(int*, int*);
typedef void (*f_GetStateVarCount)(int*, int*);
typedef void (*f_UserMod)(int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          double*,
                          double*,
                          double*,
                          double*,
                          double*,
                          const double*,
                          double*,
                          double*,
                          double*,
                          double*,
                          double**,
                          double*,
                          double*,
                          double*,
                          double*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*,
                          int*);
#endif

using SizeType = std::size_t;

SmallStrainUDSM3DLaw::SmallStrainUDSM3DLaw(const SmallStrainUDSM3DLaw& rOther)
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
    KRATOS_TRY

    for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
            mMatrixD[i][j] = rOther.mMatrixD[i][j];

    KRATOS_CATCH("")
}

ConstitutiveLaw::Pointer SmallStrainUDSM3DLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUDSM3DLaw>(*this);

    KRATOS_CATCH("")
}

SmallStrainUDSM3DLaw& SmallStrainUDSM3DLaw::operator=(SmallStrainUDSM3DLaw const& rOther)
{
    KRATOS_TRY

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

    return *this;

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::GetLawFeatures(Features& rFeatures)
{
    KRATOS_TRY

    // Set the type of law
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

    // Set strain measure required by the constitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

    // Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    KRATOS_CATCH("")
}

int SmallStrainUDSM3DLaw::Check(const Properties&   rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const ProcessInfo&  rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify Properties variables
    KRATOS_ERROR_IF(!rMaterialProperties.Has(UDSM_NAME) || rMaterialProperties[UDSM_NAME].empty())
        << "UDSM_NAME has Key zero, is not defined or has an invalid value for property"
        << rMaterialProperties.Id() << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(UDSM_NUMBER) || rMaterialProperties[UDSM_NUMBER] <= 0)
        << "UDSM_NUMBER has Key zero, is not defined or has an invalid value for property"
        << rMaterialProperties.Id() << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(IS_FORTRAN_UDSM))
        << "IS_FORTRAN_UDSM has Key zero, is not defined or has an invalid value for property"
        << rMaterialProperties.Id() << std::endl;

    return 0;
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::InitializeMaterial(const Properties&   rMaterialProperties,
                                              const GeometryType& rElementGeometry,
                                              const Vector&       rShapeFunctionsValues)
{
    KRATOS_TRY

    // loading the model
    mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    KRATOS_ERROR_IF_NOT(mIsUDSMLoaded)
        << "cannot load the specified UDSM " << rMaterialProperties[UDSM_NAME] << std::endl;

    if (rMaterialProperties[UMAT_PARAMETERS].size() != GetNumberOfMaterialParametersFromUDSM(rMaterialProperties)) {
        KRATOS_ERROR << "Number of parameters is wrong."
                     << " The UDSM gives "
                     << std::to_string(GetNumberOfMaterialParametersFromUDSM(rMaterialProperties))
                     << " while size of UMAT_PARAMETERS is "
                     << std::to_string(rMaterialProperties[UMAT_PARAMETERS].size())
                     << rMaterialProperties[UDSM_NAME] << std::endl;
    }

    ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::ResetStateVariables(const Properties& rMaterialProperties)
{
    KRATOS_TRY

    // reset state variables
    int nStateVariables = std::max(GetNumberOfStateVariablesFromUDSM(rMaterialProperties), 1);
    mStateVariables.resize(nStateVariables);
    noalias(mStateVariables) = ZeroVector(nStateVariables);

    mStateVariablesFinalized.resize(nStateVariables);
    noalias(mStateVariablesFinalized) = ZeroVector(nStateVariables);

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::ResetMaterial(const Properties&   rMaterialProperties,
                                         const GeometryType& rElementGeometry,
                                         const Vector&       rShapeFunctionsValues)
{
    KRATOS_TRY

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
    noalias(mStateVariables)          = ZeroVector(mStateVariables.size());
    noalias(mStateVariablesFinalized) = ZeroVector(mStateVariablesFinalized.size());

    mIsModelInitialized = false;

    KRATOS_CATCH(" ")
}

void SmallStrainUDSM3DLaw::SetAttributes(const Properties& rMaterialProperties)
{
    KRATOS_TRY

    if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    int IDTask = ATTRIBUTES;

    // process data
    double deltaTime = 0.0;
    double time      = 0.0;
    int    iStep     = 0;
    int    iteration = 0;

    // number of the model in the shared library (DLL)
    int modelNumber = rMaterialProperties[UDSM_NUMBER];

    // not needed:
    double bulkWater                  = 0.0;
    double excessPorePressurePrevious = 0.0;
    double excessPorePressureCurrent  = 0.0;
    double Xorigin(0.0);
    double Yorigin(0.0);
    double Zorigin(0.0);
    int    iElement          = 0;
    int    integrationNumber = 0;
    int    iPlastic          = 0;
    int    isUndr            = 0;
    int    nStateVariables   = 0;

    // variable to check if an error happend in the model:
    int                 iAbort                = 0;
    auto                nSizeProjectDirectory = static_cast<int>(mProjectDirectory.size());
    std::vector<double> StateVariablesFinalized;
    std::vector<double> StateVariables;

    const auto& MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
    pUserMod(&IDTask, &modelNumber, &isUndr, &iStep, &iteration, &iElement, &integrationNumber,
             &Xorigin, &Yorigin, &Zorigin, &time, &deltaTime, &(MaterialParameters.data()[0]),
             &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious,
             StateVariablesFinalized.data(), &(mDeltaStrainVector.data()[0]), (double**)mMatrixD,
             &bulkWater, &(mStressVector.data()[0]), &excessPorePressureCurrent,
             StateVariables.data(), &iPlastic, &nStateVariables, &mAttributes[IS_NON_SYMMETRIC],
             &mAttributes[IS_STRESS_DEPENDENT], &mAttributes[IS_TIME_DEPENDENT],
             &mAttributes[USE_TANGENT_MATRIX], mProjectDirectory.data(), &nSizeProjectDirectory, &iAbort);

    KRATOS_ERROR_IF_NOT(iAbort == 0)
        << "The specified UDSM returns an error while call UDSM with IDTASK"
        << std::to_string(IDTask) << ". UDSM" << rMaterialProperties[UDSM_NAME] << std::endl;
    KRATOS_CATCH(" ")
}

int SmallStrainUDSM3DLaw::GetNumberOfStateVariablesFromUDSM(const Properties& rMaterialProperties)
{
    KRATOS_TRY
    if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    int IDTask = NUMBER_OF_STATE_VARIABLES;

    // process data
    double deltaTime = 0.0;
    double time      = 0.0;
    int    iStep     = 0;
    int    iteration = 0;

    // number of the model in the shared library (DLL)
    int modelNumber = rMaterialProperties[UDSM_NUMBER];

    // not needed:
    double bulkWater                  = 0.0;
    double excessPorePressurePrevious = 0.0;
    double excessPorePressureCurrent  = 0.0;
    double Xorigin(0.0);
    double Yorigin(0.0);
    double Zorigin(0.0);
    int    iElement          = 0;
    int    integrationNumber = 0;
    int    iPlastic          = 0;
    int    isUndr            = 0;
    int    nStateVariables   = 0;

    // variable to check if an error happend in the model:
    int                 iAbort                = 0;
    auto                nSizeProjectDirectory = static_cast<int>(mProjectDirectory.size());
    std::vector<double> StateVariablesFinalized;
    std::vector<double> StateVariables;

    const auto& MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
    pUserMod(&IDTask, &modelNumber, &isUndr, &iStep, &iteration, &iElement, &integrationNumber,
             &Xorigin, &Yorigin, &Zorigin, &time, &deltaTime, &(MaterialParameters.data()[0]),
             &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious,
             StateVariablesFinalized.data(), &(mDeltaStrainVector.data()[0]), (double**)mMatrixD,
             &bulkWater, &(mStressVector.data()[0]), &excessPorePressureCurrent,
             StateVariables.data(), &iPlastic, &nStateVariables, &mAttributes[IS_NON_SYMMETRIC],
             &mAttributes[IS_STRESS_DEPENDENT], &mAttributes[IS_TIME_DEPENDENT],
             &mAttributes[USE_TANGENT_MATRIX], mProjectDirectory.data(), &nSizeProjectDirectory, &iAbort);

    KRATOS_ERROR_IF_NOT(iAbort == 0)
        << "The specified UDSM returns an error while call UDSM with IDTASK"
        << std::to_string(IDTask) << ". UDSM" << rMaterialProperties[UDSM_NAME] << std::endl;

    return nStateVariables;

    KRATOS_CATCH(" ")
}

SizeType SmallStrainUDSM3DLaw::GetNumberOfMaterialParametersFromUDSM(const Properties& rMaterialProperties)
{
    KRATOS_TRY

    if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    int nUDSM = rMaterialProperties[UDSM_NUMBER];
    int nParameters(0);
    pGetParamCount(&nUDSM, &nParameters);

    return static_cast<SizeType>(nParameters);

    KRATOS_CATCH("")
}

bool SmallStrainUDSM3DLaw::loadUDSM(const Properties& rMaterialProperties)
{
    KRATOS_TRY

#ifdef KRATOS_COMPILED_IN_WINDOWS
    const auto isLoaded = loadUDSMWindows(rMaterialProperties);
    return isLoaded;
#elif defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
    return loadUDSMLinux(rMaterialProperties);
#else
    KRATOS_ERROR << "loadUDSM is not supported yet for Mac OS applications" << std::endl;
#endif

    KRATOS_CATCH(" ")
}

bool SmallStrainUDSM3DLaw::loadUDSMLinux(const Properties& rMaterialProperties)
{
#ifdef KRATOS_COMPILED_IN_LINUX
    void* lib_handle;

    lib_handle = dlopen((rMaterialProperties[UDSM_NAME]).c_str(), RTLD_LAZY);
    if (!lib_handle) {
        std::string name = rMaterialProperties[UDSM_NAME];
        // check if the name of the file is based on Windows extension
        std::size_t found = name.find(".dll");
        if (found != std::string::npos) {
            // check if there is an equivalent .so file
            name.replace(found, 4, ".so");
            lib_handle = dlopen(name.c_str(), RTLD_LAZY);
        }
    }

    if (!lib_handle) {
        KRATOS_INFO("Error in loadUDSMLinux")
            << "cannot load the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
        KRATOS_ERROR << "Cannot load the specified UDSM " << rMaterialProperties[UDSM_NAME] << std::endl;
    }

    // resolve function GetParamCount address
    pGetParamCount = (f_GetParamCount)dlsym(lib_handle, "getparamcount");
    if (!pGetParamCount) {
        pGetParamCount = (f_GetParamCount)dlsym(lib_handle, "getparamcount_");
        if (!pGetParamCount) {
            KRATOS_INFO("Error in loadUDSMLinux")
                << "cannot load function GetParamCount in the specified UDSM: "
                << rMaterialProperties[UDSM_NAME] << std::endl;
            KRATOS_ERROR << "Cannot load function GetParamCount in the specified UDSM "
                         << rMaterialProperties[UDSM_NAME] << std::endl;
        }
    }

    // resolve function GetStateVarCount address
    pGetStateVarCount = (f_GetStateVarCount)dlsym(lib_handle, "getstatevarcount");

    pUserMod = (f_UserMod)dlsym(lib_handle, "user_mod");
    if (!pUserMod) {
        pUserMod = (f_UserMod)dlsym(lib_handle, "user_mod_");
        if (!pUserMod) {
            KRATOS_INFO("Error in loadUDSMLinux")
                << "cannot load function User_Mod in the specified UDSM: " << rMaterialProperties[UDSM_NAME]
                << std::endl;
            KRATOS_ERROR << "cannot load function User_Mod in the specified UDSM "
                         << rMaterialProperties[UDSM_NAME] << std::endl;
        }
    }

    return true;

#else
    KRATOS_ERROR << "loadUDSMLinux should be called in Linux applications"
                 << rMaterialProperties[UDSM_NAME] << std::endl;
#endif
}

bool SmallStrainUDSM3DLaw::loadUDSMWindows(const Properties& rMaterialProperties)
{
    KRATOS_TRY

#ifdef KRATOS_COMPILED_IN_WINDOWS

    HINSTANCE hGetProcIDDLL = LoadLibrary((rMaterialProperties[UDSM_NAME]).c_str());

    if (!hGetProcIDDLL) {
        std::string name = rMaterialProperties[UDSM_NAME];
        // check if the name of the file is based on Linux extension
        std::size_t found = name.find(".so");
        if (found != std::string::npos) {
            // check if there is an equivalent .dll file
            name.replace(found, 3, ".dll");
            hGetProcIDDLL = LoadLibrary(name.c_str());
        }
    }

    if (!hGetProcIDDLL) {
        KRATOS_INFO("Error in loadUDSMWindows")
            << "cannot load the specified UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
        KRATOS_ERROR << "cannot load the specified UDSM " << rMaterialProperties[UDSM_NAME] << std::endl;
    }

    // resolve function GetParamCount address
    pGetParamCount = (f_GetParamCount)GetProcAddress(hGetProcIDDLL, "getparamcount");
    if (!pGetParamCount) {
        // check if the dll is compiled with gfortran
        pGetParamCount = (f_GetParamCount)GetProcAddress(hGetProcIDDLL, "getparamcount_");
        if (!pGetParamCount) {
            KRATOS_INFO("Error in loadUDSMWindows")
                << "cannot load function GetParamCount in the specified UDSM: "
                << rMaterialProperties[UDSM_NAME] << std::endl;
            KRATOS_ERROR << "cannot load function GetParamCount in the specified UDSM "
                         << rMaterialProperties[UDSM_NAME] << std::endl;
        }
    }

    // resolve function GetStateVarCount address
    pGetStateVarCount = (f_GetStateVarCount)GetProcAddress(hGetProcIDDLL, "getstatevarcount");

    pUserMod = (f_UserMod)GetProcAddress(hGetProcIDDLL, "user_mod");
    if (!pUserMod) {
        // check if the dll is compiled with gfortran
        pUserMod = (f_UserMod)GetProcAddress(hGetProcIDDLL, "user_mod_");
        if (!pUserMod) {
            KRATOS_INFO("Error in loadUDSMWindows")
                << "cannot load function User_Mod in the specified UDSM: " << rMaterialProperties[UDSM_NAME]
                << std::endl;
            KRATOS_ERROR << "cannot load function User_Mod in the specified UDSM "
                         << rMaterialProperties[UDSM_NAME] << std::endl;
        }
    }
    return true;
#else
    KRATOS_ERROR << "loadUDSMWindows should be called in Windows applications"
                 << rMaterialProperties[UDSM_NAME] << std::endl;
#endif

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    CalculateMaterialResponseCauchy(rValues);
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    CalculateMaterialResponseCauchy(rValues);
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // Get Values to compute the constitutive law:
    const Flags& rOptions = rValues.GetOptions();

    // NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (rOptions.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& rStrainVector = rValues.GetStrainVector();
        CalculateCauchyGreenStrain(rValues, rStrainVector);
    }

    if (rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        // Constitutive matrix (D matrix)
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrix(rValues, rConstitutiveMatrix);
    }

    if (rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& rStressVector = rValues.GetStressVector();
        CalculateStress(rValues, rStressVector);
    }

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    const Vector& rStrainVector = rValues.GetStrainVector();

    for (unsigned int i = 0; i < mDeltaStrainVector.size(); ++i) {
        mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::SetExternalStressVector(Vector& rStressVector)
{
    KRATOS_TRY
    std::copy_n(mStressVector.begin(), rStressVector.size(), rStressVector.begin());
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::SetInternalStressVector(const Vector& rStressVector)
{
    KRATOS_TRY
    std::copy_n(rStressVector.begin(), mStressVectorFinalized.size(), mStressVectorFinalized.begin());
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY

    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer Fortran style matrix to C++ style
        for (unsigned int i = 0; i < VOIGT_SIZE_3D; i++) {
            for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[j][i];
            }
        }
    } else {
        for (unsigned int i = 0; i < VOIGT_SIZE_3D; i++) {
            for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[i][j];
            }
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CalculateConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY
    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    int IDTask = MATRIX_ELASTO_PLASTIC;
    CallUDSM(&IDTask, rValues);

    CopyConstitutiveMatrix(rValues, rConstitutiveMatrix);

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CalculateStress(ConstitutiveLaw::Parameters& rValues, Vector& rStressVector)
{
    KRATOS_TRY
    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    int IDTask = STRESS_CALCULATION;
    CallUDSM(&IDTask, rValues);

    SetExternalStressVector(rStressVector);
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::CallUDSM(int* pIDTask, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // process data
    double deltaTime = rValues.GetProcessInfo()[DELTA_TIME];
    double time      = rValues.GetProcessInfo()[TIME] - deltaTime;
    int    iStep     = rValues.GetProcessInfo()[STEP];
    int    iteration = rValues.GetProcessInfo()[NL_ITERATION_NUMBER];

    // number of the model in the shared libaray (DLL)
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    int               modelNumber         = rMaterialProperties[UDSM_NUMBER];

    // number of state variables
    auto nStateVariables = static_cast<int>(mStateVariablesFinalized.size());

    // not needed:
    double bulkWater                  = 0.0;
    double excessPorePressurePrevious = 0.0;
    double excessPorePressureCurrent  = 0.0;
    double Xorigin(0.0);
    double Yorigin(0.0);
    double Zorigin(0.0);
    int    iElement          = 0;
    int    integrationNumber = 0;
    int    iPlastic          = 0;
    int    isUndr            = 0;

    // variable to check if an error happened in the model:
    int  iAbort                = 0;
    auto nSizeProjectDirectory = static_cast<int>(mProjectDirectory.size());

    const auto& MaterialParameters = rMaterialProperties[UMAT_PARAMETERS];
    pUserMod(pIDTask, &modelNumber, &isUndr, &iStep, &iteration, &iElement, &integrationNumber,
             &Xorigin, &Yorigin, &Zorigin, &time, &deltaTime, &(MaterialParameters.data()[0]),
             &(mStressVectorFinalized.data()[0]), &excessPorePressurePrevious,
             &(mStateVariablesFinalized.data()[0]), &(mDeltaStrainVector.data()[0]),
             (double**)mMatrixD, &bulkWater, &(mStressVector.data()[0]), &excessPorePressureCurrent,
             &(mStateVariables.data()[0]), &iPlastic, &nStateVariables, &mAttributes[IS_NON_SYMMETRIC],
             &mAttributes[IS_STRESS_DEPENDENT], &mAttributes[IS_TIME_DEPENDENT],
             &mAttributes[USE_TANGENT_MATRIX], mProjectDirectory.data(), &nSizeProjectDirectory, &iAbort);

    if (iAbort != 0) {
        KRATOS_INFO("CallUDSM, iAbort !=0")
            << " iAbort: " << iAbort << " the specified UDSM returns an error while call UDSM with IDTASK: "
            << std::to_string(*pIDTask) << "."
            << " UDSM: " << rMaterialProperties[UDSM_NAME]
            << " UDSM_NUMBER: " << rMaterialProperties[UDSM_NUMBER]
            << " Parameters: " << MaterialParameters << std::endl;
        KRATOS_ERROR << "the specified UDSM returns an error while call UDSM with IDTASK: "
                     << std::to_string(*pIDTask) << ". UDSM: " << rMaterialProperties[UDSM_NAME]
                     << std::endl;
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSM3DLaw::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSM3DLaw::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSM3DLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    if (!mIsModelInitialized) {
        // stress and strain vectors must be initialized:
        const Vector& rStressVector = rValues.GetStressVector();
        SetInternalStressVector(rStressVector);

        const Vector& rStrainVector = rValues.GetStrainVector();
        SetInternalStrainVector(rStrainVector);

        int IDTask = INITIALISATION;
        CallUDSM(&IDTask, rValues);

        mIsModelInitialized = true;
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSM3DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSM3DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSM3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    UpdateInternalStrainVectorFinalized(rValues);
    mStateVariablesFinalized = mStateVariables;
    mStressVectorFinalized   = mStressVector;
}

void SmallStrainUDSM3DLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    KRATOS_TRY
    std::copy_n(rStrainVector.begin(), mStrainVectorFinalized.size(), mStrainVectorFinalized.begin());
    KRATOS_CATCH("")
}

void SmallStrainUDSM3DLaw::UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();
    this->SetInternalStrainVector(rStrainVector);
}

void SmallStrainUDSM3DLaw::CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector)
{
    const SizeType space_dimension = this->WorkingSpaceDimension();

    //-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(F.size1() != space_dimension || F.size2() != space_dimension)
        << "expected size of F " << space_dimension << "x" << space_dimension << ", got "
        << F.size1() << "x" << F.size2() << std::endl;

    Matrix E_tensor = prod(trans(F), F);
    for (unsigned int i = 0; i < space_dimension; ++i)
        E_tensor(i, i) -= 1.0;
    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

double& SmallStrainUDSM3DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                             const Variable<double>&      rThisVariable,
                                             double&                      rValue)
{
    if (rThisVariable == STRAIN_ENERGY) {
        Vector& rStrainVector = rParameterValues.GetStrainVector();
        this->CalculateCauchyGreenStrain(rParameterValues, rStrainVector);
        Vector& rStressVector = rParameterValues.GetStressVector();
        this->CalculateStress(rParameterValues, rStressVector);

        rValue = 0.5 * inner_prod(rStrainVector, rStressVector); // Strain energy = 0.5*E:C:E
    }
    return rValue;
}

Vector& SmallStrainUDSM3DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                             const Variable<Vector>&      rThisVariable,
                                             Vector&                      rValue)
{
    if (rThisVariable == STRAIN || rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {
        this->CalculateCauchyGreenStrain(rParameterValues, rValue);
    } else if (rThisVariable == STRESSES || rThisVariable == CAUCHY_STRESS_VECTOR ||
               rThisVariable == KIRCHHOFF_STRESS_VECTOR || rThisVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& rFlags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flagConstTensor = rFlags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flagStress      = rFlags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        rFlags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        rFlags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        rFlags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flagConstTensor);
        rFlags.Set(ConstitutiveLaw::COMPUTE_STRESS, flagStress);
    }

    return rValue;
}

Matrix& SmallStrainUDSM3DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                             const Variable<Matrix>&      rThisVariable,
                                             Matrix&                      rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX || rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateConstitutiveMatrix(rParameterValues, rValue);
    }
    return rValue;
}

Vector& SmallStrainUDSM3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STATE_VARIABLES) {
        if (rValue.size() != mStateVariablesFinalized.size())
            rValue.resize(mStateVariablesFinalized.size());

        noalias(rValue) = mStateVariablesFinalized;
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        if (rValue.size() != mStressVectorFinalized.size())
            rValue.resize(mStressVectorFinalized.size());

        noalias(rValue) = mStressVectorFinalized;
    }
    return rValue;
}

double& SmallStrainUDSM3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    const int index = ConstitutiveLawUtilities::GetStateVariableIndex(rThisVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "GetValue: Variable: " << rThisVariable
        << " does not exist in UDSM. Requested index: " << index << std::endl;

    rValue = mStateVariablesFinalized[index];

    return rValue;
}

int& SmallStrainUDSM3DLaw::GetValue(const Variable<int>& rThisVariable, int& rValue)
{
    if (rThisVariable == NUMBER_OF_UMAT_STATE_VARIABLES)
        rValue = static_cast<int>(mStateVariablesFinalized.size());
    return rValue;
}

void SmallStrainUDSM3DLaw::SetValue(const Variable<double>& rThisVariable,
                                    const double&           rValue,
                                    const ProcessInfo&      rCurrentProcessInfo)
{
    const int index = ConstitutiveLawUtilities::GetStateVariableIndex(rThisVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "GetValue: Variable: " << rThisVariable
        << " does not exist in UDSM. Requested index: " << index << std::endl;

    mStateVariablesFinalized[index] = rValue;
}

void SmallStrainUDSM3DLaw::SetValue(const Variable<Vector>& rThisVariable,
                                    const Vector&           rValue,
                                    const ProcessInfo&      rCurrentProcessInfo)
{
    if ((rThisVariable == STATE_VARIABLES) && (rValue.size() == mStateVariablesFinalized.size())) {
        std::copy(rValue.begin(), rValue.end(), mStateVariablesFinalized.begin());
    } else if ((rThisVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == mStressVectorFinalized.size())) {
        std::copy(rValue.begin(), rValue.end(), mStressVectorFinalized.begin());
    }
}

} // Namespace Kratos