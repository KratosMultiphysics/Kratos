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
#include <memory>
#include <string>

// Project includes
#include "includes/constitutive_law.h"
#include "includes/serializer.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{
/*
   - structure of UMATs:
   - Function to get stress, stiffness matrix
     void umat(double* STRESS,  double* STATEV,  double** DDSDDE,  double* SSE,     double* SPD,
   double* SCD, double* RPL,     double* DDSDDT,  double* DRPLDE,   double* DRPLDT,  double* STRAN,
   double* DSTRAN, double* TIME,    double* DTIME,   double* TEMP,     double* DTEMP,   double*
   PREDEF,  double* DPRED, char* MATERL,    int* NDI,        int* NSHR,        int* NTENS,      int*
   NSTATV,     double* PROPS, int* NPROPS,     double* COORDS,  double** DROT,    double* PNEWDT,
   double* CELENT,  double** DFGRD0, double** DFGRD1, int* NOEL,       int* NPT,         double*
   KSLAY,   double* KSPT,    int* KSTEP, int* KINC);

*/

class ConstitutiveLawDimension;

using pF_UMATMod = void (*)(double*       STRESS,
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

template <SizeType TVoigtSize>
class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUMATLaw : public ConstitutiveLaw
{
public:
    // The process info type definition
    using ProcessInfoType = ProcessInfo;

    // The base class ConstitutiveLaw type definition
    using BaseType = ConstitutiveLaw;

    /// The size type definition
    using SizeType = std::size_t;

    /// Pointer definition of SmallStrainUMATLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUMATLaw);

    //@}
    //@name Life Cycle
    //@{

    explicit SmallStrainUMATLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    ~SmallStrainUMATLaw() override;
    SmallStrainUMATLaw(const SmallStrainUMATLaw& rOther);
    SmallStrainUMATLaw& operator=(const SmallStrainUMATLaw& rOther);
    SmallStrainUMATLaw(SmallStrainUMATLaw&&) noexcept            = delete;
    SmallStrainUMATLaw& operator=(SmallStrainUMATLaw&&) noexcept = delete;

    /**
     * @brief Clone method
     */
    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size:
     */
    [[nodiscard]] SizeType GetStrainSize() const override;

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Infinitesimal; }

    /**
     * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
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
     * @param rVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>&      rVariable,
                           double&                      rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Vector case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Vector>&      rVariable,
                           Vector&                      rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Matrix case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Matrix>&      rVariable,
                           Matrix&                      rValue) override;

    using ConstitutiveLaw::CalculateValue;

    // @brief This function provides the place to perform checks on the completeness of the input.
    // @details It is designed to be called only once (or anyway, not often) typically at the beginning
    //          of the calculations, so to verify that nothing is missing from the input or that
    //          no common error is found.
    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
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
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void ResetMaterial(const Properties&   rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector&       rShapeFunctionsValues) override;

    using ConstitutiveLaw::GetValue;
    double& GetValue(const Variable<double>& rVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;

    using ConstitutiveLaw::SetValue;
    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    using ConstitutiveLaw::Has;
    bool Has(const Variable<Vector>& rVariable) override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    [[nodiscard]] std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected member Variables
    ///@{
    array_1d<double, TVoigtSize> mStressVector;
    array_1d<double, TVoigtSize> mStressVectorFinalized;

    array_1d<double, TVoigtSize> mDeltaStrainVector;
    array_1d<double, TVoigtSize> mStrainVectorFinalized;

    double mMatrixD[TVoigtSize][TVoigtSize];

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    SmallStrainUMATLaw();

    virtual void UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues);
    virtual void UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters& rValues);
    virtual void SetExternalStressVector(Vector& rStressVector);
    virtual void SetInternalStressVector(const Vector& rStressVector);
    virtual void SetInternalStrainVector(const Vector& rStrainVector);
    virtual void CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix);

    void CalculateConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix);
    void CalculateStress(ConstitutiveLaw::Parameters& rValues, Vector& rStressVector);

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
    pF_UMATMod mpUserMod = nullptr;

    bool mIsModelInitialized = false;
    bool mIsUMATLoaded       = false;

    std::vector<int> mProjectDirectory;

    Vector                                    mStateVariables;
    Vector                                    mStateVariablesFinalized;
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    // to load UMAT and functions
    bool loadUMAT(const Properties& rMaterialProperties);
    bool loadUMATWindows(const Properties& rMaterialProperties);
    bool loadUMATLinux(const Properties& rMaterialProperties);

    // Set number of MaterialParameters
    void CallUMAT(ConstitutiveLaw::Parameters& rValues);

    // Set state variables to the initial values
    void ResetStateVariables(const Properties& rMaterialProperties);

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SmallStrainUMATLaw

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos