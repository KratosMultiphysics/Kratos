// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_ISOTROPIC_PLASTICITY_3D_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_ISOTROPIC_PLASTICITY_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
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

/**
 * @class GenericConstitutiveLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class ConstLawIntegratorType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainIsotropicPlasticity3D
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainIsotropicPlasticity3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainIsotropicPlasticity3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainIsotropicPlasticity3D(const GenericSmallStrainIsotropicPlasticity3D &rOther)
        : ConstitutiveLaw(rOther), mPlasticDissipation(rOther.mPlasticDissipation), mThreshold(rOther.mThreshold), mPlasticStrain(rOther.mPlasticStrain), mNonConvPlasticDissipation(rOther.mNonConvPlasticDissipation), mNonConvThreshold(rOther.mNonConvThreshold), mNonConvPlasticStrain(rOther.mNonConvPlasticStrain), mUniaxialStress(rOther.mUniaxialStress)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainIsotropicPlasticity3D() override
    {
    }

    // ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override
    // {
    // }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };

    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties &rMaterialProperties,
        const GeometryType &rElementGeometry,
        const Vector &rShapeFunctionsValues,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * Finalize the material response,  called by the element in FinalizeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */
    void FinalizeMaterialResponse(ConstitutiveLaw::Parameters &rValues, const StressMeasure &rStressMeasure); //override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;
    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double> &rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector> &rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double> &rThisVariable,
        const double &rValue,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector> &rThisVariable,
        const Vector &rValue,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double &GetValue(
        const Variable<double> &rThisVariable,
        double &rValue) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector &GetValue(
        const Variable<Vector> &rThisVariable,
        Vector &rValue) override;

    /**
     * returns the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double &CalculateValue(
        ConstitutiveLaw::Parameters &rParameterValues,
        const Variable<double> &rThisVariable,
        double &rValue) override;

    Matrix &GetValue(
        const Variable<Matrix> &rThisVariable,
        Matrix &rValue) override;

    Matrix &CalculateValue(
        ConstitutiveLaw::Parameters &rParameterValues,
        const Variable<Matrix> &rThisVariable,
        Matrix &rValue) override;
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    double GetThreshold() { return mThreshold; }
    double GetPlasticDissipation() { return mPlasticDissipation; }
    Vector GetPlasticStrain() { return mPlasticStrain; }

    double GetNonConvThreshold() { return mNonConvThreshold; }
    double GetNonConvPlasticDissipation() { return mNonConvPlasticDissipation; }
    Vector GetNonConvPlasticStrain() { return mNonConvPlasticStrain; }

    void SetThreshold(const double &toThreshold) { mThreshold = toThreshold; }
    void SetPlasticDissipation(const double &toCapap) { mPlasticDissipation = toCapap; }
    void SetPlasticStrain(const Vector &Ep) { mPlasticStrain = Ep; }

    void SetNonConvThreshold(const double &toThreshold) { mNonConvThreshold = toThreshold; }
    void SetNonConvPlasticDissipation(const double &toCapap) { mNonConvPlasticDissipation = toCapap; }
    void SetNonConvPlasticStrain(const Vector &Ep) { mNonConvPlasticStrain = Ep; }


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

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

    // Converged values
    double mPlasticDissipation = 0.0;
    double mThreshold = 0.0;
    Vector mPlasticStrain = ZeroVector(6);

    // Non Converged values
    double mNonConvPlasticDissipation = 0.0;
    double mNonConvThreshold = 0.0;
    Vector mNonConvPlasticStrain = ZeroVector(6);

    // Auxiliar to print (NOTE: Alejandro do we need this now?)
    double mUniaxialStress = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues);

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor
     * @param rMaterialProperties The material properties
     */
    void CalculateElasticMatrix(
        Matrix &rElasticityTensor,
        const Properties &rMaterialProperties);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("PlasticDissipation", mPlasticDissipation);
        rSerializer.save("Threshold", mThreshold);
        rSerializer.save("PlasticStrain", mPlasticStrain);
        rSerializer.save("NonConvPlasticDissipation", mNonConvPlasticDissipation);
        rSerializer.save("NonConvThreshold", mNonConvThreshold);
        rSerializer.save("NonConvPlasticStrain", mNonConvPlasticStrain);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("PlasticDissipation", mPlasticDissipation);
        rSerializer.load("Threshold", mThreshold);
        rSerializer.load("PlasticStrain", mPlasticStrain);
        rSerializer.load("NonConvPlasticDissipation", mNonConvPlasticDissipation);
        rSerializer.load("NonConvThreshold", mNonConvThreshold);
        rSerializer.load("NonConvPlasticStrain", mNonConvPlasticStrain);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
