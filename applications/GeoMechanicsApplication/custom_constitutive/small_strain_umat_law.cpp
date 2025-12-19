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

#include "custom_constitutive/small_strain_umat_law.h"
#include "constitutive_law_dimension.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

#ifdef KRATOS_COMPILED_IN_WINDOWS
#include "windows.hpp"
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
#include <dlfcn.h>
#endif

namespace Kratos
{
using namespace std::string_literals;

#ifdef KRATOS_COMPILED_IN_WINDOWS
using f_UMATMod = void(__stdcall*)(double*       STRESS,
                                   double*       STATEV,
                                   double**      DDSDDE,
                                   double*       SSE,
                                   double*       SPD,
                                   double*       SCD,
                                   double*       RPL,
                                   double*       DDSDDT,
                                   double*       DRPLDE,
                                   double*       DRPLDT,
                                   double*       STRAN,
                                   double*       DSTRAN,
                                   double*       TIME,
                                   double*       DTIME,
                                   double*       TEMP,
                                   double*       DTEMP,
                                   double*       PREDEF,
                                   double*       DPRED,
                                   char*         MATERL,
                                   int*          NDI,
                                   int*          NSHR,
                                   int*          NTENS,
                                   int*          NSTATV,
                                   const double* PROPS,
                                   int*          NPROPS,
                                   double*       COORDS,
                                   double**      DROT,
                                   double*       PNEWDT,
                                   double*       CELENT,
                                   double**      DFGRD0,
                                   double**      DFGRD1,
                                   int*          NOEL,
                                   int*          NPT,
                                   double*       KSLAY,
                                   double*       KSPT,
                                   int*          KSTEP,
                                   int*          KINC);
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
using f_UMATMod = void (*)(double*       STRESS,
                           double*       STATEV,
                           double**      DDSDDE,
                           double*       SSE,
                           double*       SPD,
                           double*       SCD,
                           double*       rpl,
                           double*       ddsddt,
                           double*       drplde,
                           double*       drpldt,
                           double*       stran,
                           double*       dstran,
                           double*       time,
                           double*       dtime,
                           double*       temp,
                           double*       dtemp,
                           double*       predef,
                           double*       dpred,
                           char*         materl,
                           int*          ndi,
                           int*          nshr,
                           int*          ntens,
                           int*          nstatv,
                           const double* props,
                           int*          nprops,
                           double*       coords,
                           double**      drot,
                           double*       pnewdt,
                           double*       celent,
                           double**      dfgrd0,
                           double**      dfgrd1,
                           int*          noel,
                           int*          npt,
                           double*       kslay,
                           double*       kspt,
                           int*          kstep,
                           int*          kinc);
#endif

template <SizeType TVoigtSize>
SmallStrainUMATLaw<TVoigtSize>::SmallStrainUMATLaw() = default;

template <SizeType TVoigtSize>
SmallStrainUMATLaw<TVoigtSize>::SmallStrainUMATLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : ConstitutiveLaw{}, mpConstitutiveDimension(std::move(pConstitutiveDimension))
{
}

template <SizeType TVoigtSize>
SmallStrainUMATLaw<TVoigtSize>::~SmallStrainUMATLaw() = default;

template <SizeType TVoigtSize>
SmallStrainUMATLaw<TVoigtSize>::SmallStrainUMATLaw(const SmallStrainUMATLaw& rOther)
    : ConstitutiveLaw(rOther),
      mStressVector(rOther.mStressVector),
      mStressVectorFinalized(rOther.mStressVectorFinalized),
      mDeltaStrainVector(rOther.mDeltaStrainVector),
      mStrainVectorFinalized(rOther.mStrainVectorFinalized),
      mIsModelInitialized(rOther.mIsModelInitialized),
      mIsUMATLoaded(rOther.mIsUMATLoaded),
      mStateVariables(rOther.mStateVariables),
      mStateVariablesFinalized(rOther.mStateVariablesFinalized),
      mpConstitutiveDimension(rOther.mpConstitutiveDimension->Clone())
{
    KRATOS_TRY

    for (unsigned int i = 0; i < TVoigtSize; ++i)
        for (unsigned int j = 0; j < TVoigtSize; ++j)
            mMatrixD[i][j] = rOther.mMatrixD[i][j];

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
ConstitutiveLaw::Pointer SmallStrainUMATLaw<TVoigtSize>::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUMATLaw>(*this);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
SmallStrainUMATLaw<TVoigtSize>& SmallStrainUMATLaw<TVoigtSize>::operator=(const SmallStrainUMATLaw<TVoigtSize>& rOther)
{
    KRATOS_TRY

    ConstitutiveLaw::operator=(rOther);
    this->mIsModelInitialized      = rOther.mIsModelInitialized;
    this->mIsUMATLoaded            = rOther.mIsUMATLoaded;
    this->mStateVariables          = rOther.mStateVariables;
    this->mStateVariablesFinalized = rOther.mStateVariablesFinalized;
    this->mStressVector            = rOther.mStressVector;
    this->mStressVectorFinalized   = rOther.mStressVectorFinalized;
    this->mDeltaStrainVector       = rOther.mDeltaStrainVector;
    this->mStrainVectorFinalized   = rOther.mStrainVectorFinalized;
    this->mpConstitutiveDimension  = rOther.mpConstitutiveDimension->Clone();

    for (unsigned int i = 0; i < TVoigtSize; ++i)
        for (unsigned int j = 0; j < TVoigtSize; ++j)
            this->mMatrixD[i][j] = rOther.mMatrixD[i][j];

    return *this;

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set(mpConstitutiveDimension->GetSpatialType());
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);

    rFeatures.mOptions.Set(ISOTROPIC);

    // Set strain measure required by the constitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

    // Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();
}

template <SizeType TVoigtSize>
int SmallStrainUMATLaw<TVoigtSize>::Check(const Properties&   rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const ProcessInfo&  rCurrentProcessInfo) const
{
    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllExclusive);
    check_properties.CheckAvailabilityAndNotEmpty(UDSM_NAME);
    check_properties.CheckAvailability(IS_FORTRAN_UDSM);

    return 0;
}

template <SizeType TVoigtSize>
SizeType SmallStrainUMATLaw<TVoigtSize>::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

template <SizeType TVoigtSize>
SizeType SmallStrainUMATLaw<TVoigtSize>::GetStrainSize() const
{
    // In other constitutive laws, we use mpConstitutiveDimension->GetStrainSize() here, but
    // due to the C/Fortran interface, we need the VoigtSize to be known compile time.
    // Therefore, we return the template argument TVoigtSize here.
    return TVoigtSize;
}

template <SizeType TVoigtSize>
ConstitutiveLaw::StrainMeasure SmallStrainUMATLaw<TVoigtSize>::GetStrainMeasure()
{
    return StrainMeasure_Infinitesimal;
}

template <SizeType TVoigtSize>
ConstitutiveLaw::StressMeasure SmallStrainUMATLaw<TVoigtSize>::GetStressMeasure()
{
    return StressMeasure_Cauchy;
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::InitializeMaterial(const Properties&   rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const Vector&       rShapeFunctionsValues)

{
    KRATOS_TRY
    // we need to check if the model is loaded or not
    mIsUMATLoaded = loadUMAT(rMaterialProperties);

    if (!mIsUMATLoaded)
        KRATOS_ERROR << "cannot load the specified UMAT" << rMaterialProperties[UDSM_NAME] << std::endl;

    ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::ResetStateVariables(const Properties& rMaterialProperties)
{
    KRATOS_TRY
    // reset state variables

    const auto& state_variables   = rMaterialProperties[STATE_VARIABLES];
    const auto  n_state_variables = state_variables.size();

    mStateVariables.resize(n_state_variables);
    mStateVariablesFinalized.resize(n_state_variables);

    noalias(mStateVariables)          = state_variables;
    noalias(mStateVariablesFinalized) = state_variables;

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("InitializedModel", mIsModelInitialized);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("StateVariablesFinalized", mStateVariablesFinalized);
    rSerializer.save("ConstitutitiveLawDimension", mpConstitutiveDimension);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("InitializedModel", mIsModelInitialized);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("StateVariablesFinalized", mStateVariablesFinalized);
    rSerializer.load("ConstitutitiveLawDimension", mpConstitutiveDimension);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::ResetMaterial(const Properties&   rMaterialProperties,
                                                   const GeometryType& rElementGeometry,
                                                   const Vector&       rShapeFunctionsValues)
{
    KRATOS_TRY

    // reset state variables
    ResetStateVariables(rMaterialProperties);

    // set stress vectors:
    noalias(mStressVector)          = ZeroVector(mStressVector.size());
    noalias(mStressVectorFinalized) = ZeroVector(mStressVectorFinalized.size());

    // set strain vectors:
    noalias(mDeltaStrainVector)     = ZeroVector(mDeltaStrainVector.size());
    noalias(mStrainVectorFinalized) = ZeroVector(mStrainVectorFinalized.size());

    for (unsigned int i = 0; i < TVoigtSize; ++i)
        for (unsigned int j = 0; j < TVoigtSize; ++j)
            mMatrixD[i][j] = 0.0;

    mIsModelInitialized = false;

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
bool SmallStrainUMATLaw<TVoigtSize>::loadUMAT(const Properties& rMaterialProperties)
{
    KRATOS_TRY

#ifdef KRATOS_COMPILED_IN_WINDOWS
    return loadUMATWindows(rMaterialProperties);
#elif defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
    return loadUMATLinux(rMaterialProperties);
#else
    KRATOS_ERROR << "loadUMAT is not supported yet for Mac OS applications" << std::endl;
#endif

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
bool SmallStrainUMATLaw<TVoigtSize>::loadUMATLinux(const Properties& rMaterialProperties)
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
        KRATOS_ERROR << "cannot load the specified UMAT " << rMaterialProperties[UDSM_NAME] << std::endl;
        return false;
    }
    if (rMaterialProperties[IS_FORTRAN_UDSM]) {
        mpUserMod = (f_UMATMod)dlsym(lib_handle, "umat_");
    } else {
        mpUserMod = (f_UMATMod)dlsym(lib_handle, "umat");
    }

    if (!mpUserMod) {
        KRATOS_ERROR << "cannot load function User_Mod in the specified UMAT "
                     << rMaterialProperties[UDSM_NAME] << std::endl;
        return false;
    }

    return true;

#else
    KRATOS_ERROR << "loadUMATLinux should be called in Linux applications"
                 << rMaterialProperties[UDSM_NAME] << std::endl;
#endif
}

template <SizeType TVoigtSize>
bool SmallStrainUMATLaw<TVoigtSize>::loadUMATWindows(const Properties& rMaterialProperties)
{
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
        KRATOS_INFO("Error in loadUMATWindows")
            << "cannot load the specified UMAT: " << rMaterialProperties[UDSM_NAME] << std::endl;
        KRATOS_ERROR << "cannot load the specified UMAT " << rMaterialProperties[UDSM_NAME] << std::endl;
    }

    mpUserMod = (f_UMATMod)GetProcAddress(hGetProcIDDLL, "umat");
    if (!mpUserMod) {
        KRATOS_INFO("Error in loadUMATWindows")
            << "cannot load function umat in the specified UMAT: " << rMaterialProperties[UDSM_NAME]
            << std::endl;
        KRATOS_ERROR << "cannot load function umat in the specified UMAT "
                     << rMaterialProperties[UDSM_NAME] << std::endl;
    }

    return true;
#else
    KRATOS_ERROR << "loadUMATWindows should be called in Windows applications"
                 << rMaterialProperties[UDSM_NAME] << std::endl;
    return false;
#endif
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // Get Values to compute the constitutive law:
    const Flags& r_options = rValues.GetOptions();

    KRATOS_DEBUG_ERROR_IF(r_options.IsDefined(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN) &&
                          r_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
        << "The SmallStrainUMATLaw needs an element provided strain" << std::endl;

    KRATOS_ERROR_IF(!rValues.IsSetStrainVector() || rValues.GetStrainVector().size() != GetStrainSize())
        << "Constitutive laws in the geomechanics application need a valid provided strain" << std::endl;

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculateStress(rValues, r_stress_vector);
    }

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        // Constitutive matrix (D matrix)
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrix(rValues, r_constitutive_matrix);
    }

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& r_strain_vector = rValues.GetStrainVector();

    for (unsigned int i = 0; i < TVoigtSize; ++i) {
        mDeltaStrainVector[i] = r_strain_vector(i) - mStrainVectorFinalized[i];
    }
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::SetExternalStressVector(Vector& rStressVector)
{
    std::copy_n(mStressVector.begin(), TVoigtSize, rStressVector.begin());
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::SetInternalStressVector(const Vector& rStressVector)
{
    std::copy_n(rStressVector.begin(), TVoigtSize, mStressVectorFinalized.begin());
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::copy_n(rStrainVector.begin(), TVoigtSize, mStrainVectorFinalized.begin());
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues,
                                                            Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < TVoigtSize; i++) {
            for (unsigned int j = 0; j < TVoigtSize; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[j][i];
            }
        }
    } else {
        for (unsigned int i = 0; i < TVoigtSize; i++) {
            for (unsigned int j = 0; j < TVoigtSize; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[i][j];
            }
        }
    }
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CalculateConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues,
                                                                 Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY

    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    CallUMAT(rValues);

    CopyConstitutiveMatrix(rValues, rConstitutiveMatrix);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CalculateStress(ConstitutiveLaw::Parameters& rValues, Vector& rStressVector)
{
    KRATOS_TRY

    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    CallUMAT(rValues);

    SetExternalStressVector(rStressVector);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::CallUMAT(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // process data
    double deltaTime = rValues.GetProcessInfo()[DELTA_TIME];
    double time      = rValues.GetProcessInfo()[TIME] - deltaTime;
    int    iStep     = rValues.GetProcessInfo()[STEP];
    int    iteration = rValues.GetProcessInfo()[NL_ITERATION_NUMBER];

    // number of the model in the shared libaray (DLL)

    // number of state variables
    auto nStateVariables = static_cast<int>(mStateVariablesFinalized.size());

    // not needed:
    int    iElement          = 0;
    int    integrationNumber = 0;
    double SPD; // specific plastic dissipation
    double SSE; // ?
    double SCD; // ?
    char   materialName;

    int ndi   = mpConstitutiveDimension->GetNumberOfNormalComponents();
    int ntens = TVoigtSize;
    int nshr  = ntens - ndi;

    // stresses and state variables in the beginning of the steps needs to be given:
    mStressVector   = mStressVectorFinalized;
    mStateVariables = mStateVariablesFinalized;

    // variable to check if an error happened in the model:
    const auto& MaterialParameters = rValues.GetMaterialProperties()[UMAT_PARAMETERS];
    auto        nProperties        = static_cast<int>(MaterialParameters.size());
    mpUserMod(&(mStressVector.data()[0]), &(mStateVariables.data()[0]), (double**)mMatrixD, &SSE,
              &SPD, &SCD, nullptr, nullptr, nullptr, nullptr, &(mStrainVectorFinalized.data()[0]),
              &(mDeltaStrainVector.data()[0]), &time, &deltaTime, nullptr, nullptr, nullptr,
              nullptr, &materialName, &ndi, &nshr, &ntens, &nStateVariables,
              &(MaterialParameters.data()[0]), &nProperties, nullptr, nullptr, nullptr, nullptr,
              nullptr, nullptr, &iElement, &integrationNumber, nullptr, nullptr, &iStep, &iteration);

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    if (!mIsModelInitialized) {
        // stress and strain vectors must be initialized:
        const Vector& r_stress_vector = rValues.GetStressVector();
        const Vector& r_strain_vector = rValues.GetStrainVector();

        SetInternalStressVector(r_stress_vector);

        SetInternalStrainVector(r_strain_vector);

        CallUMAT(rValues);
        mIsModelInitialized = true;
    }

    KRATOS_CATCH("")
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    UpdateInternalStrainVectorFinalized(rValues);
    mStateVariablesFinalized = mStateVariables;
    mStressVectorFinalized   = mStressVector;
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();
    this->SetInternalStrainVector(rStrainVector);
}

template <SizeType TVoigtSize>
double& SmallStrainUMATLaw<TVoigtSize>::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                       const Variable<double>& rVariable,
                                                       double&                 rValue)
{
    if (rVariable == STRAIN_ENERGY) {
        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector&       r_stress_vector = rParameterValues.GetStressVector();
        this->CalculateStress(rParameterValues, r_stress_vector);

        rValue = 0.5 * inner_prod(r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return rValue;
}

template <SizeType TVoigtSize>
Vector& SmallStrainUMATLaw<TVoigtSize>::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                       const Variable<Vector>& rVariable,
                                                       Vector&                 rValue)
{
    if (rVariable == STRESSES || rVariable == CAUCHY_STRESS_VECTOR ||
        rVariable == KIRCHHOFF_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& rFlags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flagConstTensor = rFlags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flagStress      = rFlags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        rFlags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        rFlags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        SmallStrainUMATLaw::CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        rFlags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flagConstTensor);
        rFlags.Set(ConstitutiveLaw::COMPUTE_STRESS, flagStress);
    }

    return rValue;
}

template <SizeType TVoigtSize>
Matrix& SmallStrainUMATLaw<TVoigtSize>::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                       const Variable<Matrix>& rVariable,
                                                       Matrix&                 rValue)
{
    if (rVariable == CONSTITUTIVE_MATRIX || rVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateConstitutiveMatrix(rParameterValues, rValue);
    }

    return rValue;
}

template <SizeType TVoigtSize>
Vector& SmallStrainUMATLaw<TVoigtSize>::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STATE_VARIABLES) {
        rValue.resize(mStateVariablesFinalized.size());
        noalias(rValue) = mStateVariablesFinalized;
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue.resize(mStressVectorFinalized.size());
        noalias(rValue) = mStressVectorFinalized;
    }

    return rValue;
}

template <SizeType TVoigtSize>
double& SmallStrainUMATLaw<TVoigtSize>::GetValue(const Variable<double>& rVariable, double& rValue)
{
    int index = ConstitutiveLawUtilities::GetStateVariableIndex(rVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "GetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

    rValue = mStateVariablesFinalized[index];

    return rValue;
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::SetValue(const Variable<double>& rVariable,
                                              const double&           rValue,
                                              const ProcessInfo&      rCurrentProcessInfo)
{
    const int index = ConstitutiveLawUtilities::GetStateVariableIndex(rVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "SetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

    mStateVariablesFinalized[index] = rValue;
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::SetValue(const Variable<Vector>& rVariable,
                                              const Vector&           rValue,
                                              const ProcessInfo&      rCurrentProcessInfo)
{
    if ((rVariable == STATE_VARIABLES) && (rValue.size() == mStateVariablesFinalized.size())) {
        std::ranges::copy(rValue, mStateVariablesFinalized.begin());
    } else if ((rVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == TVoigtSize)) {
        std::ranges::copy_n(rValue.begin(), TVoigtSize, mStressVectorFinalized.begin());
    }
}

template <SizeType TVoigtSize>
bool SmallStrainUMATLaw<TVoigtSize>::Has(const Variable<Vector>& rVariable)
{
    return rVariable == STATE_VARIABLES || rVariable == CAUCHY_STRESS_VECTOR;
}

template <SizeType TVoigtSize>
[[nodiscard]] std::string SmallStrainUMATLaw<TVoigtSize>::Info() const
{
    return "SmallStrainUMATLaw"s;
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

template <SizeType TVoigtSize>
void SmallStrainUMATLaw<TVoigtSize>::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUMATLaw Data";
}

template class SmallStrainUMATLaw<VOIGT_SIZE_3D>;
template class SmallStrainUMATLaw<VOIGT_SIZE_2D_INTERFACE>;
template class SmallStrainUMATLaw<VOIGT_SIZE_3D_INTERFACE>;

} // Namespace Kratos
