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

#if !defined(KRATOS_VISCOUS_GENERALIZED_MAXWELL_H_INCLUDED)
#define KRATOS_VISCOUS_GENERALIZED_MAXWELL_H_INCLUDED

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
 * @class ViscousGeneralizedMaxwell3D
 * @ingroup StructuralMechanicsApplication
 * @brief This is a viscous law using Maxwell formulation
 * @details The definition of a maxwell material can be found in https://en.wikipedia.org/wiki/Maxwell_material
 * This definition consists in a spring and a damper in serial
 *
 *           -----^^^^^^-----------[------
 *                Spring (K)   Damper (C)
 *
 * @author Alejandro Cornejo & Lucia Barbu
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ViscousGeneralizedMaxwell3D
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    typedef ConstitutiveLaw BaseType;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ViscousGeneralizedMaxwell3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    ViscousGeneralizedMaxwell3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ViscousGeneralizedMaxwell3D>(*this);
    }

    /**
    * Copy constructor.
    */
    ViscousGeneralizedMaxwell3D(const ViscousGeneralizedMaxwell3D &rOther)
        : BaseType(rOther),
          mPrevStressVector(rOther.mPrevStressVector),
          mPrevStrainVector(rOther.mPrevStrainVector),
          mNonConvPrevStressVector(rOther.mNonConvPrevStressVector),
          mNonConvPrevStrainVector(rOther.mNonConvPrevStrainVector)
    {
    }

    /**
    * Destructor.
    */
    ~ViscousGeneralizedMaxwell3D() override
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
        const Properties &rMaterialProperties,
        const GeometryType &rElementGeometry,
        const Vector &rShapeFunctionsValues,
        const ProcessInfo &rCurrentProcessInfo) override;

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
    Vector mPrevStrainVector = ZeroVector(6);

    // Non Converged values
    Vector mNonConvPrevStressVector = ZeroVector(6);
    Vector mNonConvPrevStrainVector = ZeroVector(6);

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

    Vector GetPreviousStrainVector() { return mPrevStrainVector; }
    void SetPreviousStrainVector(Vector toStrain) { mPrevStrainVector = toStrain; }
    Vector GetNonConvPreviousStrainVector() { return mNonConvPrevStrainVector; }
    void SetNonConvPreviousStrainVector(Vector toStrain) { mNonConvPrevStrainVector = toStrain; }

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor
     * @param rMaterialProperties The material properties
     */
    void CalculateElasticMatrix(
        Matrix& rElasticityTensor,
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
        rSerializer.save("PrevStressVector", mPrevStressVector);
        rSerializer.save("PrevStrainVector", mPrevStrainVector);
        rSerializer.save("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.save("NonConvPrevStrainVector", mNonConvPrevStrainVector);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("PrevStressVector", mPrevStressVector);
        rSerializer.load("PrevStrainVector", mPrevStrainVector);
        rSerializer.load("NonConvPrevStressVector", mNonConvPrevStressVector);
        rSerializer.load("NonConvPrevStrainVector", mNonConvPrevStrainVector);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
