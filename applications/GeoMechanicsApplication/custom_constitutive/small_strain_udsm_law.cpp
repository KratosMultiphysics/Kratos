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

#include <algorithm>
#include <type_traits>

#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/small_strain_udsm_law.hpp"
#include "custom_utilities/constitutive_law_utilities.h"

#ifdef KRATOS_COMPILED_IN_WINDOWS
#include "windows.hpp"
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
#include <dlfcn.h>
#endif

namespace
{

using namespace Kratos;

// See section 16.2 "Implementation of User Defined (UD) soil Models in calculations program"
// of the Plaxis documentation for the array size of `Props`
constexpr auto props_size = SizeType{50};

array_1d<double, props_size> MakePropsVector(const Vector& rUMatParameters)
{
    KRATOS_DEBUG_ERROR_IF(rUMatParameters.size() > props_size)
        << "Number of UMAT_PARAMETERS (" << rUMatParameters.size()
        << ") exceeds the maximum number of " << props_size << "\n";

    auto result = array_1d<double, props_size>{props_size, 0.0};
    std::copy(rUMatParameters.begin(), rUMatParameters.end(), result.begin());
    return result;
}

} // namespace

namespace Kratos
{

// calling convention (__cdecl, __stdcall, ...)
// __stdcall is the convention used by the WinAPI
#ifdef KRATOS_COMPILED_IN_WINDOWS
using f_GetParamCount    = void(__stdcall*)(int*, int*);
using f_GetStateVarCount = void(__stdcall*)(int*, int*);
using f_UserMod          = void(__stdcall*)(int*,
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
using f_GetParamCount    = void (*)(int*, int*);
using f_GetStateVarCount = void (*)(int*, int*);
using f_UserMod          = void (*)(int*,
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

SmallStrainUDSMLaw::~SmallStrainUDSMLaw() = default;

SmallStrainUDSMLaw::SmallStrainUDSMLaw(std::unique_ptr<ConstitutiveLawDimension> pDimension)
    : mpDimension(std::move(pDimension))
{
    ResetConstitutiveMatrix();
}

SmallStrainUDSMLaw::SmallStrainUDSMLaw(SmallStrainUDSMLaw&&) noexcept            = default;
SmallStrainUDSMLaw& SmallStrainUDSMLaw::operator=(SmallStrainUDSMLaw&&) noexcept = default;

ConstitutiveLaw::Pointer SmallStrainUDSMLaw::Clone() const
{
    auto pResult = std::make_shared<SmallStrainUDSMLaw>();
    CloneDataMembersTo(*pResult);
    return pResult;
}

void SmallStrainUDSMLaw::GetLawFeatures(Features& rFeatures)
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

SizeType SmallStrainUDSMLaw::WorkingSpaceDimension() { return mpDimension->GetDimension(); }

SizeType SmallStrainUDSMLaw::GetStrainSize() const { return mpDimension->GetStrainSize(); }

ConstitutiveLaw::StrainMeasure SmallStrainUDSMLaw::GetStrainMeasure()
{
    return StrainMeasure_Infinitesimal;
}

ConstitutiveLaw::StressMeasure SmallStrainUDSMLaw::GetStressMeasure()
{
    return StressMeasure_Cauchy;
}

int SmallStrainUDSMLaw::Check(const Properties&   rMaterialProperties,
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

void SmallStrainUDSMLaw::InitializeMaterial(const Properties&   rMaterialProperties,
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

void SmallStrainUDSMLaw::ResetStateVariables(const Properties& rMaterialProperties)
{
    KRATOS_TRY

    // reset state variables
    int n_state_variables = std::max(GetNumberOfStateVariablesFromUDSM(rMaterialProperties), 1);
    mStateVariables.resize(n_state_variables);
    noalias(mStateVariables) = ZeroVector(n_state_variables);

    mStateVariablesFinalized.resize(n_state_variables);
    noalias(mStateVariablesFinalized) = ZeroVector(n_state_variables);

    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::ResetConstitutiveMatrix()
{
    for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i) {
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j) {
            mMatrixD[i][j] = 0.0;
        }
    }
}

void SmallStrainUDSMLaw::ResetMaterial(const Properties&   rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector&       rShapeFunctionsValues)
{
    KRATOS_TRY

    // reset state variables
    SetAttributes(rMaterialProperties);
    ResetStateVariables(rMaterialProperties);

    // set stress vectors:
    noalias(mStressVector) = ZeroVector(mStressVector.size());
    mSig0.clear();

    // set strain vectors:
    mDeltaStrainVector.clear();
    noalias(mStrainVectorFinalized) = ZeroVector(mStrainVectorFinalized.size());

    ResetConstitutiveMatrix();

    // state variables
    noalias(mStateVariables)          = ZeroVector(mStateVariables.size());
    noalias(mStateVariablesFinalized) = ZeroVector(mStateVariablesFinalized.size());

    mIsModelInitialized = false;

    KRATOS_CATCH(" ")
}

void SmallStrainUDSMLaw::SetAttributes(const Properties& rMaterialProperties)
{
    KRATOS_TRY

    if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    auto task_id_as_int = static_cast<int>(UdsmTaskId::ATTRIBUTES);

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

    const auto umat_parameters = MakePropsVector(rMaterialProperties[UMAT_PARAMETERS]);

    mpUserMod(&task_id_as_int, &modelNumber, &isUndr, &iStep, &iteration, &iElement, &integrationNumber,
              &Xorigin, &Yorigin, &Zorigin, &time, &deltaTime, &(umat_parameters.data()[0]),
              &(mSig0.data()[0]), &excessPorePressurePrevious, StateVariablesFinalized.data(),
              &(mDeltaStrainVector.data()[0]), (double**)mMatrixD, &bulkWater,
              &(mStressVector.data()[0]), &excessPorePressureCurrent, StateVariables.data(), &iPlastic,
              &nStateVariables, &mAttributes[IS_NON_SYMMETRIC], &mAttributes[IS_STRESS_DEPENDENT],
              &mAttributes[IS_TIME_DEPENDENT], &mAttributes[USE_TANGENT_MATRIX],
              mProjectDirectory.data(), &nSizeProjectDirectory, &iAbort);

    KRATOS_ERROR_IF_NOT(iAbort == 0)
        << "The specified UDSM returns an error while call UDSM with IDTASK" << task_id_as_int
        << ". UDSM" << rMaterialProperties[UDSM_NAME] << std::endl;
    KRATOS_CATCH(" ")
}

int SmallStrainUDSMLaw::GetNumberOfStateVariablesFromUDSM(const Properties& rMaterialProperties)
{
    KRATOS_TRY
    if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    auto task_id_as_int = static_cast<int>(UdsmTaskId::NUMBER_OF_STATE_VARIABLES);

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

    // variable to check if an error occurred in the model:
    int                 iAbort                = 0;
    auto                nSizeProjectDirectory = static_cast<int>(mProjectDirectory.size());
    std::vector<double> StateVariablesFinalized;
    std::vector<double> StateVariables;

    const auto umat_parameters = MakePropsVector(rMaterialProperties[UMAT_PARAMETERS]);

    mpUserMod(&task_id_as_int, &modelNumber, &isUndr, &iStep, &iteration, &iElement, &integrationNumber,
              &Xorigin, &Yorigin, &Zorigin, &time, &deltaTime, &(umat_parameters.data()[0]),
              &(mSig0.data()[0]), &excessPorePressurePrevious, StateVariablesFinalized.data(),
              &(mDeltaStrainVector.data()[0]), (double**)mMatrixD, &bulkWater,
              &(mStressVector.data()[0]), &excessPorePressureCurrent, StateVariables.data(), &iPlastic,
              &nStateVariables, &mAttributes[IS_NON_SYMMETRIC], &mAttributes[IS_STRESS_DEPENDENT],
              &mAttributes[IS_TIME_DEPENDENT], &mAttributes[USE_TANGENT_MATRIX],
              mProjectDirectory.data(), &nSizeProjectDirectory, &iAbort);

    KRATOS_ERROR_IF_NOT(iAbort == 0)
        << "The specified UDSM returns an error while call UDSM with IDTASK" << task_id_as_int
        << ". UDSM" << rMaterialProperties[UDSM_NAME] << std::endl;

    return nStateVariables;

    KRATOS_CATCH(" ")
}

SizeType SmallStrainUDSMLaw::GetNumberOfMaterialParametersFromUDSM(const Properties& rMaterialProperties)
{
    KRATOS_TRY

    if (!mIsUDSMLoaded) mIsUDSMLoaded = loadUDSM(rMaterialProperties);

    int nUDSM = rMaterialProperties[UDSM_NUMBER];
    int nParameters(0);
    mpGetParamCount(&nUDSM, &nParameters);

    return static_cast<SizeType>(nParameters);

    KRATOS_CATCH("")
}

bool SmallStrainUDSMLaw::loadUDSM(const Properties& rMaterialProperties)
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

bool SmallStrainUDSMLaw::loadUDSMLinux(const Properties& rMaterialProperties)
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
    mpGetParamCount = (f_GetParamCount)dlsym(lib_handle, "getparamcount");
    if (!mpGetParamCount) {
        mpGetParamCount = (f_GetParamCount)dlsym(lib_handle, "getparamcount_");
        if (!mpGetParamCount) {
            KRATOS_INFO("Error in loadUDSMLinux")
                << "cannot load function GetParamCount in the specified UDSM: "
                << rMaterialProperties[UDSM_NAME] << std::endl;
            KRATOS_ERROR << "Cannot load function GetParamCount in the specified UDSM "
                         << rMaterialProperties[UDSM_NAME] << std::endl;
        }
    }

    // resolve function GetStateVarCount address
    mpGetStateVarCount = (f_GetStateVarCount)dlsym(lib_handle, "getstatevarcount");

    mpUserMod = (f_UserMod)dlsym(lib_handle, "user_mod");
    if (!mpUserMod) {
        mpUserMod = (f_UserMod)dlsym(lib_handle, "user_mod_");
        if (!mpUserMod) {
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

bool SmallStrainUDSMLaw::loadUDSMWindows(const Properties& rMaterialProperties)
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
    mpGetParamCount = (f_GetParamCount)GetProcAddress(hGetProcIDDLL, "getparamcount");
    if (!mpGetParamCount) {
        // check if the dll is compiled with gfortran
        mpGetParamCount = (f_GetParamCount)GetProcAddress(hGetProcIDDLL, "getparamcount_");
        if (!mpGetParamCount) {
            KRATOS_INFO("Error in loadUDSMWindows")
                << "cannot load function GetParamCount in the specified UDSM: "
                << rMaterialProperties[UDSM_NAME] << std::endl;
            KRATOS_ERROR << "cannot load function GetParamCount in the specified UDSM "
                         << rMaterialProperties[UDSM_NAME] << std::endl;
        }
    }

    // resolve function GetStateVarCount address
    mpGetStateVarCount = (f_GetStateVarCount)GetProcAddress(hGetProcIDDLL, "getstatevarcount");

    mpUserMod = (f_UserMod)GetProcAddress(hGetProcIDDLL, "user_mod");
    if (!mpUserMod) {
        // check if the dll is compiled with gfortran
        mpUserMod = (f_UserMod)GetProcAddress(hGetProcIDDLL, "user_mod_");
        if (!mpUserMod) {
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

void SmallStrainUDSMLaw::CalculateMaterialResponsePK1(Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    CalculateMaterialResponseCauchy(rValues);
    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    KRATOS_TRY
    CalculateMaterialResponseCauchy(rValues);
    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_TRY

    // Get Values to compute the constitutive law:
    const Flags& rOptions = rValues.GetOptions();

    KRATOS_DEBUG_ERROR_IF(rOptions.IsDefined(USE_ELEMENT_PROVIDED_STRAIN) && rOptions.IsNot(USE_ELEMENT_PROVIDED_STRAIN))
        << "The SmallStrainUDSMLaw needs an element provided strain" << std::endl;

    KRATOS_ERROR_IF(!rValues.IsSetStrainVector() || rValues.GetStrainVector().size() != GetStrainSize())
        << "Constitutive laws in the geomechanics application need a valid provided strain" << std::endl;

    if (rOptions.Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
        CalculateConstitutiveMatrix(rValues, rValues.GetConstitutiveMatrix());
    }

    if (rOptions.Is(COMPUTE_STRESS)) {
        CalculateStress(rValues, rValues.GetStressVector());
    }

    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::UpdateInternalDeltaStrainVector(Parameters& rValues)
{
    const auto& r_strain_vector = rValues.GetStrainVector();
    std::transform(r_strain_vector.begin(),
                   r_strain_vector.begin() + static_cast<Vector::difference_type>(GetStrainSize()),
                   mStrainVectorFinalized.begin(), mDeltaStrainVector.begin(), std::minus<>{});
}

void SmallStrainUDSMLaw::SetExternalStressVector(Vector& rStressVector)
{
    std::copy_n(mStressVector.begin(), GetStrainSize(), rStressVector.begin());
}

void SmallStrainUDSMLaw::SetInternalStressVector(const Vector& rStressVector)
{
    std::copy(rStressVector.begin(), rStressVector.end(), mSig0.begin());
}

void SmallStrainUDSMLaw::CopyConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY

    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer Fortran style matrix to C++ style
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[j][i];
            }
        }
    } else {
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[i][j];
            }
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::CloneDataMembersTo(SmallStrainUDSMLaw& rDestination) const
{
    rDestination.mStressVector          = mStressVector;
    rDestination.mDeltaStrainVector     = mDeltaStrainVector;
    rDestination.mStrainVectorFinalized = mStrainVectorFinalized;
    for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i) {
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j) {
            rDestination.mMatrixD[i][j] = mMatrixD[i][j];
        }
    }
    rDestination.mpGetParamCount          = mpGetParamCount;
    rDestination.mpGetStateVarCount       = mpGetStateVarCount;
    rDestination.mpUserMod                = mpUserMod;
    rDestination.mIsModelInitialized      = mIsModelInitialized;
    rDestination.mIsUDSMLoaded            = mIsUDSMLoaded;
    rDestination.mAttributes              = mAttributes;
    rDestination.mProjectDirectory        = mProjectDirectory;
    rDestination.mStateVariables          = mStateVariables;
    rDestination.mStateVariablesFinalized = mStateVariablesFinalized;
    rDestination.mSig0                    = mSig0;
    rDestination.mpDimension              = mpDimension ? mpDimension->Clone() : nullptr;
}

void SmallStrainUDSMLaw::CalculateConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY
    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    CallUDSM(UdsmTaskId::MATRIX_ELASTO_PLASTIC, rValues);

    CopyConstitutiveMatrix(rValues, rConstitutiveMatrix);

    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::CalculateStress(Parameters& rValues, Vector& rStressVector)
{
    KRATOS_TRY
    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    CallUDSM(UdsmTaskId::STRESS_CALCULATION, rValues);

    SetExternalStressVector(rStressVector);
    KRATOS_CATCH("")
}

array_1d<double, SmallStrainUDSMLaw::Sig0Size>& SmallStrainUDSMLaw::GetSig0() { return mSig0; }

void SmallStrainUDSMLaw::CallUDSM(UdsmTaskId TaskId, Parameters& rValues)
{
    KRATOS_TRY

    auto task_id_as_int = static_cast<int>(TaskId);

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

    const auto umat_parameters = MakePropsVector(rMaterialProperties[UMAT_PARAMETERS]);

    mpUserMod(&task_id_as_int, &modelNumber, &isUndr, &iStep, &iteration, &iElement, &integrationNumber,
              &Xorigin, &Yorigin, &Zorigin, &time, &deltaTime, &(umat_parameters.data()[0]),
              &(mSig0.data()[0]), &excessPorePressurePrevious, &(mStateVariablesFinalized.data()[0]),
              &(mDeltaStrainVector.data()[0]), (double**)mMatrixD, &bulkWater,
              &(mStressVector.data()[0]), &excessPorePressureCurrent, &(mStateVariables.data()[0]),
              &iPlastic, &nStateVariables, &mAttributes[IS_NON_SYMMETRIC],
              &mAttributes[IS_STRESS_DEPENDENT], &mAttributes[IS_TIME_DEPENDENT],
              &mAttributes[USE_TANGENT_MATRIX], mProjectDirectory.data(), &nSizeProjectDirectory, &iAbort);

    if (iAbort != 0) {
        KRATOS_INFO("CallUDSM, iAbort !=0")
            << " iAbort: " << iAbort
            << " the specified UDSM returns an error while call UDSM with IDTASK: " << task_id_as_int << "."
            << " UDSM: " << rMaterialProperties[UDSM_NAME]
            << " UDSM_NUMBER: " << rMaterialProperties[UDSM_NUMBER]
            << " Parameters: " << rMaterialProperties[UMAT_PARAMETERS] << std::endl;
        KRATOS_ERROR << "the specified UDSM returns an error while call UDSM with IDTASK: " << task_id_as_int
                     << ". UDSM: " << rMaterialProperties[UDSM_NAME] << std::endl;
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::InitializeMaterialResponsePK1(Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSMLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSMLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSMLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_TRY

    if (!mIsModelInitialized) {
        SetInternalStressVector(rValues.GetStressVector());
        SetInternalStrainVector(rValues.GetStrainVector());

        CallUDSM(UdsmTaskId::INITIALISATION, rValues);

        mIsModelInitialized = true;
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSMLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSMLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSMLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

void SmallStrainUDSMLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    UpdateInternalStrainVectorFinalized(rValues);
    mStateVariablesFinalized = mStateVariables;
    std::copy(mStressVector.begin(), mStressVector.end(), mSig0.begin());
}

void SmallStrainUDSMLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::copy_n(rStrainVector.begin(), GetStrainSize(), mStrainVectorFinalized.begin());
}

void SmallStrainUDSMLaw::UpdateInternalStrainVectorFinalized(Parameters& rValues)
{
    SetInternalStrainVector(rValues.GetStrainVector());
}

double& SmallStrainUDSMLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rVariable, double& rValue)
{
    if (rVariable == STRAIN_ENERGY) {
        CalculateStress(rParameterValues, rParameterValues.GetStressVector());
        rValue = 0.5 * inner_prod(rParameterValues.GetStrainVector(), rParameterValues.GetStressVector()); // Strain energy = 0.5*E:C:E
    }
    return rValue;
}

Vector& SmallStrainUDSMLaw::CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STRESSES || rVariable == CAUCHY_STRESS_VECTOR ||
        rVariable == KIRCHHOFF_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& rFlags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flagConstTensor = rFlags.Is(COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flagStress      = rFlags.Is(COMPUTE_STRESS);

        rFlags.Set(COMPUTE_CONSTITUTIVE_TENSOR, true);
        rFlags.Set(COMPUTE_STRESS, true);

        // We compute the stress
        CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        rFlags.Set(COMPUTE_CONSTITUTIVE_TENSOR, flagConstTensor);
        rFlags.Set(COMPUTE_STRESS, flagStress);
    }

    return rValue;
}

Matrix& SmallStrainUDSMLaw::CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rVariable, Matrix& rValue)
{
    if (rVariable == CONSTITUTIVE_MATRIX || rVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        CalculateConstitutiveMatrix(rParameterValues, rValue);
    }
    return rValue;
}

Vector& SmallStrainUDSMLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STATE_VARIABLES) {
        rValue.resize(mStateVariablesFinalized.size());
        noalias(rValue) = mStateVariablesFinalized;
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue.resize(GetStrainSize());
        std::copy_n(mSig0.begin(), GetStrainSize(), rValue.begin());
    }
    return rValue;
}

double& SmallStrainUDSMLaw::GetValue(const Variable<double>& rVariable, double& rValue)
{
    const int index = ConstitutiveLawUtilities::GetStateVariableIndex(rVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "GetValue: Variable: " << rVariable
        << " does not exist in UDSM. Requested index: " << index << std::endl;

    rValue = mStateVariablesFinalized[index];

    return rValue;
}

void SmallStrainUDSMLaw::SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    const int index = ConstitutiveLawUtilities::GetStateVariableIndex(rVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "GetValue: Variable: " << rVariable
        << " does not exist in UDSM. Requested index: " << index << std::endl;

    mStateVariablesFinalized[index] = rValue;
}

void SmallStrainUDSMLaw::SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ((rVariable == STATE_VARIABLES) && (rValue.size() == mStateVariablesFinalized.size())) {
        std::copy(rValue.begin(), rValue.end(), mStateVariablesFinalized.begin());
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        KRATOS_ERROR_IF(rValue.size() != GetStrainSize())
            << "Failed to set stress vector: expected one with " << GetStrainSize()
            << " components, but got one with " << rValue.size() << "components\n";
        std::copy(rValue.begin(), rValue.end(), mSig0.begin());
    }
}

std::string SmallStrainUDSMLaw::Info() const { return "SmallStrainUDSMLaw"; }

void SmallStrainUDSMLaw::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

void SmallStrainUDSMLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUDSMLaw Data";
}

void SmallStrainUDSMLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)

    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    // Saving `mMatrixD` was missing and I wasn't able to find an easy way to add it, since it is a
    // C-style array. When the type of `mMatrixD` is modified to a proper matrix type, this will
    // become a no-brainer.

    // Also, it doesn't make sense to save function pointers pointing into the shared library, since
    // it may no longer be loaded when we restore this constitutive law. Therefore, member `load`
    // keeps the initial null values. By setting the 'is UDSM loaded' flag to false we enforce the
    // shared library to be loaded again before using it.
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
    rSerializer.save("IsUDSMLoaded", false);
    rSerializer.save("Attributes", mAttributes);
    rSerializer.save("ProjectDirectory", mProjectDirectory);
    rSerializer.save("StateVariables", mStateVariables);
    rSerializer.save("StateVariablesFinalized", mStateVariablesFinalized);
    rSerializer.save("Sig0", mSig0);
    rSerializer.save("Dimension", mpDimension);
}

void SmallStrainUDSMLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)

    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("DeltaStrainVector", mDeltaStrainVector);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    // Loading `mMatrixD` was missing and cannot be added yet, since saving it is rather
    // complicated (see also the comment in member `save`)

    // Also the function pointers cannot be restored. They will keep their initial null values.

    rSerializer.load("IsModelInitialized", mIsModelInitialized);
    rSerializer.load("IsUDSMLoaded", mIsUDSMLoaded);
    rSerializer.load("Attributes", mAttributes);
    rSerializer.load("ProjectDirectory", mProjectDirectory);
    rSerializer.load("StateVariables", mStateVariables);
    rSerializer.load("StateVariablesFinalized", mStateVariablesFinalized);
    rSerializer.load("Sig0", mSig0);
    rSerializer.load("Dimension", mpDimension);
}

// Instances of this class cannot be copied, but they can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<SmallStrainUDSMLaw>);
static_assert(!std::is_copy_assignable_v<SmallStrainUDSMLaw>);
static_assert(std::is_move_constructible_v<SmallStrainUDSMLaw>);
static_assert(std::is_move_assignable_v<SmallStrainUDSMLaw>);

} // Namespace Kratos