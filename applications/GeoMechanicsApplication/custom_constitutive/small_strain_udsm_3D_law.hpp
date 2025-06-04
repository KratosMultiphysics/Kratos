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

#pragma once

#include <iostream>
#include <string>

#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

class ConstitutiveLawDimension;

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

using pF_GetParamCount    = void (*)(int*, int*);
using pF_GetStateVarCount = void (*)(int*, int*);
using pF_UserMod          = void (*)(int*,
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

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSMLaw : public ConstitutiveLaw
{
public:
    using SizeType = std::size_t;

    // See section 16.2 "Implementation of User Defined (UD) soil Models in calculations program"
    // of the Plaxis documentation for the array sizes
    static constexpr SizeType Sig0Size                  = 20;
    static constexpr SizeType StressVectorSize          = 6;
    static constexpr SizeType StrainIncrementVectorSize = 12;

    /// Pointer definition of SmallStrainUDSMLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSMLaw);

    explicit SmallStrainUDSMLaw(std::unique_ptr<ConstitutiveLawDimension> pDimension = nullptr);
    ~SmallStrainUDSMLaw() override;

    // This constitutive law cannot be copied. Use member function `Clone` instead.
    SmallStrainUDSMLaw(const SmallStrainUDSMLaw&)            = delete;
    SmallStrainUDSMLaw& operator=(const SmallStrainUDSMLaw&) = delete;

    // This constitutive law can be moved
    SmallStrainUDSMLaw(SmallStrainUDSMLaw&&) noexcept;
    SmallStrainUDSMLaw& operator=(SmallStrainUDSMLaw&&) noexcept;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    void GetLawFeatures(Features& rFeatures) override;

    SizeType WorkingSpaceDimension() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    StrainMeasure GetStrainMeasure() override;
    StressMeasure GetStressMeasure() override;

    void CalculateMaterialResponsePK1(Parameters& rValues) override;
    void CalculateMaterialResponsePK2(Parameters& rValues) override;
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;
    void FinalizeMaterialResponsePK1(Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(Parameters& rValues) override;
    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override;

    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rVariable, double& rValue) override;
    Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rVariable, Vector& rValue) override;
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rVariable, Matrix& rValue) override;
    using ConstitutiveLaw::CalculateValue;

    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;

    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    void InitializeMaterialResponseCauchy(Parameters& rValues) override;
    void InitializeMaterialResponsePK1(Parameters& rValues) override;
    void InitializeMaterialResponsePK2(Parameters& rValues) override;
    void InitializeMaterialResponseKirchhoff(Parameters& rValues) override;

    void ResetMaterial(const Properties&   rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector&       rShapeFunctionsValues) override;

    double& GetValue(const Variable<double>& rVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    using ConstitutiveLaw::GetValue;

    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using ConstitutiveLaw::SetValue;

    [[nodiscard]] std::string Info() const override;
    void                      PrintInfo(std::ostream& rOStream) const override;
    void                      PrintData(std::ostream& rOStream) const override;

protected:
    enum IDTASK : int {
        INITIALISATION = 1,
        STRESS_CALCULATION,
        MATRIX_ELASTO_PLASTIC,
        NUMBER_OF_STATE_VARIABLES,
        ATTRIBUTES,
        MATRIX_ELASTIC
    };

    enum ATTRIBUTE : int {
        IS_NON_SYMMETRIC,
        IS_STRESS_DEPENDENT,
        IS_TIME_DEPENDENT,
        USE_TANGENT_MATRIX
    };

    array_1d<double, VOIGT_SIZE_3D> mStressVector{VOIGT_SIZE_3D, 0.0};

    array_1d<double, StrainIncrementVectorSize> mDeltaStrainVector{StrainIncrementVectorSize, 0.0};
    array_1d<double, VOIGT_SIZE_3D>             mStrainVectorFinalized{VOIGT_SIZE_3D, 0.0};

    double mMatrixD[VOIGT_SIZE_3D][VOIGT_SIZE_3D];

    virtual void UpdateInternalDeltaStrainVector(Parameters& rValues);
    virtual void UpdateInternalStrainVectorFinalized(Parameters& rValues);
    virtual void SetExternalStressVector(Vector& rStressVector);
    virtual void SetInternalStressVector(const Vector& rStressVector);
    virtual void SetInternalStrainVector(const Vector& rStrainVector);
    virtual void CopyConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix);

    void CloneDataMembersTo(SmallStrainUDSMLaw& rDestination) const;

    void CalculateConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix);
    void CalculateStress(Parameters& rValues, Vector& rStressVector);

    // returns 1 if the stiffness matrix of the material is non-symmetric
    int getIsNonSymmetric();

    // returns 1 if the stiffness matrix of the material is stress dependent
    int getIsStressDependent();

    // returns 1 if material is time dependent
    int getIsTimeDependent();

    // returns 1 if the stiffness matrix of the material is tangential
    int getUseTangentMatrix();

    array_1d<double, Sig0Size>& GetSig0();

private:
    pF_GetParamCount    mpGetParamCount    = nullptr;
    pF_GetStateVarCount mpGetStateVarCount = nullptr;
    pF_UserMod          mpUserMod          = nullptr;

    bool mIsModelInitialized = false;
    bool mIsUDSMLoaded       = false;

    array_1d<int, 4> mAttributes;

    std::vector<int> mProjectDirectory;

    Vector mStateVariables;
    Vector mStateVariablesFinalized;

    array_1d<double, Sig0Size> mSig0{Sig0Size, 0.0};

    std::unique_ptr<ConstitutiveLawDimension> mpDimension;

    // to load UDSM and functions
    bool loadUDSM(const Properties& rMaterialProperties);
    bool loadUDSMWindows(const Properties& rMaterialProperties);
    bool loadUDSMLinux(const Properties& rMaterialProperties);

    void CallUDSM(int* IDTask, ConstitutiveLaw::Parameters& rValues);

    void ResetStateVariables(const Properties& rMaterialProperties);

    void SetAttributes(const Properties& rMaterialProperties);

    int GetNumberOfStateVariablesFromUDSM(const Properties& rMaterialProperties);

    SizeType          GetNumberOfMaterialParametersFromUDSM(const Properties& rMaterialProperties);
    [[nodiscard]] int GetStateVariableIndex(const Variable<double>& rVariable) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos