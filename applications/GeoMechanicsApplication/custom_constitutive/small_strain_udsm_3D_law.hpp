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

// System includes
#include "includes/define.h"
#include <iostream>
#include <string>

// Project includes
#include "includes/constitutive_law.h"
#include "includes/serializer.h"

// Application includes
#include "geo_mechanics_application_variables.h"

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

typedef void (*pF_GetParamCount)(int*, int*);
typedef void (*pF_GetStateVarCount)(int*, int*);
typedef void (*pF_UserMod)(int*,
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

///@addtogroup ConstitutiveModelsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM3DLaw : public ConstitutiveLaw
{
public:
    // The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw BaseType;

    /// The size type definition
    typedef std::size_t SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_3D;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = VOIGT_SIZE_3D;

    /// Pointer definition of SmallStrainUDSM3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM3DLaw);

    //@}
    //@name Life Cycle
    //@{

    //----------------------------------------------------------------------------------------
    /**
     * @brief Default constructor.
     */
    SmallStrainUDSM3DLaw() = default;

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    SmallStrainUDSM3DLaw(SmallStrainUDSM3DLaw const& rOther);

    /**
     * @brief Destructor.
     */
    ~SmallStrainUDSM3DLaw() override = default;

    // Assignment operator:
    SmallStrainUDSM3DLaw& operator=(SmallStrainUDSM3DLaw const& rOther);

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override { return Dimension; }

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override { return VoigtSize; }

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Infinitesimal; }

    /**
     * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in Voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override { return StressMeasure_Cauchy; }

    /**
     * @brief Computes the material response:
     * @details PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response:
     * @details Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response:
     * @details Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Updates the material response:
     * @details Cauchy stresses and Internal Variables
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>&      rThisVariable,
                           double&                      rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Vector case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Vector>&      rThisVariable,
                           Vector&                      rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Matrix case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Matrix>&      rThisVariable,
                           Matrix&                      rValue) override;

    using ConstitutiveLaw::CalculateValue;

    // @brief This function provides the place to perform checks on the completeness of the input.
    // @details It is designed to be called only once (or anyway, not often) typically at the beginning
    //          of the calculations, so to verify that nothing is missing from the input or that
    //          no common error is found.
    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    /**
     * Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the strain vector
     * @param rValues The internal values of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    virtual void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector);

    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void ResetMaterial(const Properties&   rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector&       rShapeFunctionsValues) override;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    int&    GetValue(const Variable<int>& rThisVariable, int& rValue) override;
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override { return "SmallStrainUDSM3DLaw"; }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "SmallStrainUDSM3DLaw Data";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
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

    ///@}
    ///@name Protected member Variables
    ///@{
    array_1d<double, VOIGT_SIZE_3D> mStressVector;
    array_1d<double, VOIGT_SIZE_3D> mStressVectorFinalized;

    array_1d<double, VOIGT_SIZE_3D> mDeltaStrainVector;
    array_1d<double, VOIGT_SIZE_3D> mStrainVectorFinalized;

    double mMatrixD[VOIGT_SIZE_3D][VOIGT_SIZE_3D];

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{
    virtual void UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues);
    virtual void UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters& rValues);
    virtual void SetExternalStressVector(Vector& rStressVector);
    virtual void SetInternalStressVector(const Vector& rStressVector);
    virtual void SetInternalStrainVector(const Vector& rStrainVector);
    virtual void CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix);

    void CalculateConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix);
    void CalculateStress(ConstitutiveLaw::Parameters& rValues, Vector& rStressVector);

    // returns 1 if the stiffness matrix of the material is non-symmetric
    int getIsNonSymmetric() { return mAttributes[IS_NON_SYMMETRIC]; }

    // returns 1 if the stiffness matrix of the material is stress dependent
    int getIsStressDependent() { return mAttributes[IS_STRESS_DEPENDENT]; }

    // returns 1 if material is time dependent
    int getIsTimeDependent() { return mAttributes[IS_TIME_DEPENDENT]; }

    // returns 1 if the stiffness matrix of the material is tangential
    int getUseTangentMatrix() { return mAttributes[USE_TANGENT_MATRIX]; }

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    pF_GetParamCount    pGetParamCount;
    pF_GetStateVarCount pGetStateVarCount;
    pF_UserMod          pUserMod;

    bool mIsModelInitialized = false;
    bool mIsUDSMLoaded       = false;

    // AttributesUDSM mAttributes;
    array_1d<int, 4> mAttributes;

    std::vector<int> mProjectDirectory;

    Vector mStateVariables;
    Vector mStateVariablesFinalized;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    // to load UDSM and functions
    bool loadUDSM(const Properties& rMaterialProperties);
    bool loadUDSMWindows(const Properties& rMaterialProperties);
    bool loadUDSMLinux(const Properties& rMaterialProperties);

    // Set number of MaterialParameters
    void CallUDSM(int* IDTask, ConstitutiveLaw::Parameters& rValues);

    // Set state variables to zero
    void ResetStateVariables(const Properties& rMaterialProperties);

    // set number of StateVariables
    void SetAttributes(const Properties& rMaterialProperties);

    // set number of StateVariables
    int GetNumberOfStateVariablesFromUDSM(const Properties& rMaterialProperties);

    // get number of MaterialParameters
    SizeType GetNumberOfMaterialParametersFromUDSM(const Properties& rMaterialProperties);
    int      GetStateVariableIndex(const Variable<double>& rThisVariable) const;

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("InitializedModel", mIsModelInitialized);
        rSerializer.save("Attributes", mAttributes);
        rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
        rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
        rSerializer.save("StateVariablesFinalized", mStateVariablesFinalized);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("InitializedModel", mIsModelInitialized);
        rSerializer.load("Attributes", mAttributes);
        rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
        rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
        rSerializer.load("StateVariablesFinalized", mStateVariablesFinalized);
    }

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SmallStrainUDSM3DLaw

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos