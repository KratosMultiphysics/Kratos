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
#include "custom_utilities/constitutive_law_utilities.hpp"

#ifdef KRATOS_COMPILED_IN_WINDOWS
#include "windows.hpp"
#endif

#ifdef KRATOS_COMPILED_IN_LINUX
#include <dlfcn.h>
#endif

namespace Kratos
{

#ifdef KRATOS_COMPILED_IN_WINDOWS
typedef void(__stdcall* f_UMATMod)(double*       STRESS,
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
typedef void (*f_UMATMod)(double*       STRESS,
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

SmallStrainUMAT3DLaw::SmallStrainUMAT3DLaw(const SmallStrainUMAT3DLaw& rOther)
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
    KRATOS_TRY

    for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
            mMatrixD[i][j] = rOther.mMatrixD[i][j];

    KRATOS_CATCH("")
}

ConstitutiveLaw::Pointer SmallStrainUMAT3DLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUMAT3DLaw>(*this);

    KRATOS_CATCH("")
}

SmallStrainUMAT3DLaw& SmallStrainUMAT3DLaw::operator=(const SmallStrainUMAT3DLaw& rOther)
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

    for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
            this->mMatrixD[i][j] = rOther.mMatrixD[i][j];

    return *this;

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);

    rFeatures.mOptions.Set(ISOTROPIC);

    // Set strain measure required by the constitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

    // Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();
}

int SmallStrainUMAT3DLaw::Check(const Properties&   rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const ProcessInfo&  rCurrentProcessInfo) const
{
    // Verify Properties variables
    if (!rMaterialProperties.Has(UDSM_NAME) || rMaterialProperties[UDSM_NAME] == "")
        KRATOS_ERROR << "UDSM_NAME has Key zero, is not defined for property"
                     << rMaterialProperties.Id() << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(IS_FORTRAN_UDSM))
        << "IS_FORTRAN_UDSM is not defined for property" << rMaterialProperties.Id() << std::endl;

    return 0;
}

void SmallStrainUMAT3DLaw::InitializeMaterial(const Properties&   rMaterialProperties,
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

void SmallStrainUMAT3DLaw::ResetStateVariables(const Properties& rMaterialProperties)
{
    KRATOS_TRY
    // reset state variables

    const auto& StateVariables  = rMaterialProperties[STATE_VARIABLES];
    const auto  nStateVariables = StateVariables.size();

    mStateVariables.resize(nStateVariables);
    mStateVariablesFinalized.resize(nStateVariables);

    noalias(mStateVariables)          = StateVariables;
    noalias(mStateVariablesFinalized) = StateVariables;

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::ResetMaterial(const Properties&   rMaterialProperties,
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

    for (unsigned int i = 0; i < VOIGT_SIZE_3D; ++i)
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; ++j)
            mMatrixD[i][j] = 0.0;

    mIsModelInitialized = false;

    KRATOS_CATCH("")
}

bool SmallStrainUMAT3DLaw::loadUMAT(const Properties& rMaterialProperties)
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

bool SmallStrainUMAT3DLaw::loadUMATLinux(const Properties& rMaterialProperties)
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
        pUserMod = (f_UMATMod)dlsym(lib_handle, "umat_");
    } else {
        pUserMod = (f_UMATMod)dlsym(lib_handle, "umat");
    }

    if (!pUserMod) {
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

bool SmallStrainUMAT3DLaw::loadUMATWindows(const Properties& rMaterialProperties)
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

    pUserMod = (f_UMATMod)GetProcAddress(hGetProcIDDLL, "umat");
    if (!pUserMod) {
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

void SmallStrainUMAT3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // Get Values to compute the constitutive law:
    const Flags& rOptions = rValues.GetOptions();

    // NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (rOptions.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& rStrainVector = rValues.GetStrainVector();
        CalculateCauchyGreenStrain(rValues, rStrainVector);
    }

    if (rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& rStressVector = rValues.GetStressVector();
        CalculateStress(rValues, rStressVector);
    }

    if (rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        // Constitutive matrix (D matrix)
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrix(rValues, rConstitutiveMatrix);
    }

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();

    for (unsigned int i = 0; i < mDeltaStrainVector.size(); ++i) {
        mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
    }
}

void SmallStrainUMAT3DLaw::SetExternalStressVector(Vector& rStressVector)
{
    for (unsigned int i = 0; i < rStressVector.size(); ++i) {
        rStressVector(i) = mStressVector[i];
    }
}

void SmallStrainUMAT3DLaw::SetInternalStressVector(const Vector& rStressVector)
{
    for (unsigned int i = 0; i < mStressVectorFinalized.size(); ++i) {
        mStressVectorFinalized[i] = rStressVector(i);
    }
}

void SmallStrainUMAT3DLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    for (unsigned int i = 0; i < mStrainVectorFinalized.size(); ++i) {
        mStrainVectorFinalized[i] = rStrainVector(i);
    }
}

void SmallStrainUMAT3DLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
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
}

void SmallStrainUMAT3DLaw::CalculateConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY

    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    CallUMAT(rValues);

    CopyConstitutiveMatrix(rValues, rConstitutiveMatrix);

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::CalculateStress(ConstitutiveLaw::Parameters& rValues, Vector& rStressVector)
{
    KRATOS_TRY

    // update strain vector
    UpdateInternalDeltaStrainVector(rValues);

    CallUMAT(rValues);

    SetExternalStressVector(rStressVector);

    KRATOS_CATCH("")
}

void SmallStrainUMAT3DLaw::CallUMAT(ConstitutiveLaw::Parameters& rValues)
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

    int ndi   = N_DIM_3D;
    int nshr  = 3;
    int ntens = VOIGT_SIZE_3D;

    // stresses and state variables in the beginning of the steps needs to be given:
    mStressVector   = mStressVectorFinalized;
    mStateVariables = mStateVariablesFinalized;

    // variable to check if an error happened in the model:
    const auto& MaterialParameters = rValues.GetMaterialProperties()[UMAT_PARAMETERS];
    auto        nProperties        = static_cast<int>(MaterialParameters.size());
    pUserMod(&(mStressVector.data()[0]), &(mStateVariables.data()[0]), (double**)mMatrixD, &SSE,
             &SPD, &SCD, nullptr, nullptr, nullptr, nullptr, &(mStrainVectorFinalized.data()[0]),
             &(mDeltaStrainVector.data()[0]), &time, &deltaTime, nullptr, nullptr, nullptr, nullptr,
             &materialName, &ndi, &nshr, &ntens, &nStateVariables, &(MaterialParameters.data()[0]),
             &nProperties, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &iElement,
             &integrationNumber, nullptr, nullptr, &iStep, &iteration);

    KRATOS_CATCH("")
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
    KRATOS_TRY

    if (!mIsModelInitialized) {
        // stress and strain vectors must be initialized:
        const Vector& rStressVector = rValues.GetStressVector();
        const Vector& rStrainVector = rValues.GetStrainVector();

        SetInternalStressVector(rStressVector);

        SetInternalStrainVector(rStrainVector);

        CallUMAT(rValues);
        mIsModelInitialized = true;
    }

    KRATOS_CATCH("")
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

void SmallStrainUMAT3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    UpdateInternalStrainVectorFinalized(rValues);
    mStateVariablesFinalized = mStateVariables;
    mStressVectorFinalized   = mStressVector;
}

void SmallStrainUMAT3DLaw::UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();
    this->SetInternalStrainVector(rStrainVector);
}

void SmallStrainUMAT3DLaw::CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector)
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

double& SmallStrainUMAT3DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                             const Variable<double>&      rThisVariable,
                                             double&                      rValue)
{
    Vector& rStrainVector = rParameterValues.GetStrainVector();
    Vector& rStressVector = rParameterValues.GetStressVector();

    if (rThisVariable == STRAIN_ENERGY) {
        this->CalculateCauchyGreenStrain(rParameterValues, rStrainVector);
        this->CalculateStress(rParameterValues, rStressVector);

        rValue = 0.5 * inner_prod(rStrainVector, rStressVector); // Strain energy = 0.5*E:C:E
    }

    return rValue;
}

Vector& SmallStrainUMAT3DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
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
        SmallStrainUMAT3DLaw::CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        rFlags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flagConstTensor);
        rFlags.Set(ConstitutiveLaw::COMPUTE_STRESS, flagStress);
    }

    return rValue;
}

Matrix& SmallStrainUMAT3DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                             const Variable<Matrix>&      rThisVariable,
                                             Matrix&                      rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX || rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateConstitutiveMatrix(rParameterValues, rValue);
    }

    return rValue;
}

Vector& SmallStrainUMAT3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
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

double& SmallStrainUMAT3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    int index = ConstitutiveLawUtilities::GetStateVariableIndex(rThisVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "GetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

    rValue = mStateVariablesFinalized[index];

    return rValue;
}

int& SmallStrainUMAT3DLaw::GetValue(const Variable<int>& rThisVariable, int& rValue)
{
    if (rThisVariable == NUMBER_OF_UMAT_STATE_VARIABLES)
        rValue = static_cast<int>(mStateVariablesFinalized.size());

    return rValue;
}

void SmallStrainUMAT3DLaw::SetValue(const Variable<double>& rThisVariable,
                                    const double&           rValue,
                                    const ProcessInfo&      rCurrentProcessInfo)
{
    const int index = ConstitutiveLawUtilities::GetStateVariableIndex(rThisVariable);

    KRATOS_DEBUG_ERROR_IF(index < 0 || index > (static_cast<int>(mStateVariablesFinalized.size()) - 1))
        << "SetValue: State variable does not exist in UDSM. Requested index: " << index << std::endl;

    mStateVariablesFinalized[index] = rValue;
}

void SmallStrainUMAT3DLaw::SetValue(const Variable<Vector>& rThisVariable,
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
