// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo&  Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_VISCOUS_GENERALIZED_KELVIN_H_INCLUDED)
#define KRATOS_VISCOUS_GENERALIZED_KELVIN_H_INCLUDED

// System includes

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
 * @class ViscousGeneralizedKelvin3D
 * @ingroup StructuralMechanicsApplication
 * @brief This is a constitutive law that reproduces the behaviour of viscous Kelvin material
 * @details The definition of the Kelvin-Voigt material can be founf in Wikipedia https://en.wikipedia.org/wiki/Kelvin%E2%80%93Voigt_material
 * The Kelvin–Voigt model, also called the Voigt model, can be represented by a purely viscous damper and purely elastic spring
 *
 *                   Spring K
 *                  |--^^^^--|
 *             -----          ------
 *                  |---[----|
 *                   Damper C
 *
 * @author Alejandro Cornejo&  Lucia Barbu
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ViscousGeneralizedKelvin3D
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base class
    typedef ConstitutiveLaw BaseType;

    /// The index definition
    typedef std::size_t IndexType;

    /// The size definition
    typedef std::size_t SizeType;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ViscousGeneralizedKelvin3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    ViscousGeneralizedKelvin3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ViscousGeneralizedKelvin3D>(*this);
    }

    /**
    * Copy constructor.
    */
    ViscousGeneralizedKelvin3D(const ViscousGeneralizedKelvin3D& rOther)
        : ConstitutiveLaw(rOther)
    {
    }
    /**
    * Destructor.
    */
    ~ViscousGeneralizedKelvin3D() override
    {
    }

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
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

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
    Vector mPrevStressVector = ZeroVector(6);
    Vector mPrevInelasticStrainVector = ZeroVector(6);

    // Non Converged values
    Vector mNonConvPrevStressVector = ZeroVector(6);
    Vector mNonConvPrevInelasticStrainVector = ZeroVector(6);

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    Vector GetPreviousStressVector() { return mPrevStressVector; }
    void SetPreviousStressVector(Vector toStress) { mPrevStressVector = toStress; }
    Vector GetNonConvPreviousStressVector() { return mNonConvPrevStressVector; }
    void SetNonConvPreviousStressVector(Vector toStress) { mNonConvPrevStressVector = toStress; }

    Vector GetPreviousInelasticStrainVector() { return mPrevInelasticStrainVector; }
    void SetPreviousInelasticStrainVector(Vector toStrain) { mPrevInelasticStrainVector = toStrain; }
    Vector GetNonConvPreviousInelasticStrainVector() { return mNonConvPrevInelasticStrainVector; }
    void SetNonConvPreviousInelasticStrainVector(Vector toStrain) { mNonConvPrevInelasticStrainVector = toStrain; }

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor
     * @param rMaterialProperties The material properties
     */
    void CalculateElasticMatrix(
        Matrix& rElasticityTensor,
        const Properties& rMaterialProperties);

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("PrevStressVector", mPrevStressVector);
        rSerializer.save("PrevInelasticStrainVector", mPrevInelasticStrainVector);
        rSerializer.save("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.save("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("PrevStressVector", mPrevStressVector);
        rSerializer.load("PrevInelasticStrainVector", mPrevInelasticStrainVector);
        rSerializer.load("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.load("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
