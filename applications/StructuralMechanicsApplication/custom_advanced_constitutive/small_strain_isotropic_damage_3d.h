// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//  Collaborator:

#if !defined(KRATOS_SMALL_STRAIN_ISOTROPIC_DAMAGE_3D_LAW_H_INCLUDED)
#define KRATOS_SMALL_STRAIN_ISOTROPIC_DAMAGE_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "elastic_isotropic_3d.h"

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
 * @class SmallStrainIsotropicDamage3D
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a damage with hardening/softening constitutive law in 3D
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * - INFINITY_YIELD_STRESS
 * - HARDENING_MODULI_VECTOR: List of hardening modules (only two branches considered)
 * @warning Valid for small strains
 * @note
 * @author Marcelo Raschi
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallStrainIsotropicDamage3D
    : public ElasticIsotropic3D

{
public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    // Counted pointer of LinearIsotropicDamage3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainIsotropicDamage3D);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainIsotropicDamage3D();

    /**
     * @brief Default constructor.
     */
    SmallStrainIsotropicDamage3D(const SmallStrainIsotropicDamage3D& rOther);

    /**
     * @brief Default constructor.
     */
    ~SmallStrainIsotropicDamage3D() override;

    /**
     * @brief Clone function
     * @return A pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    ///@}
    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (bool)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (bool)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    bool& GetValue(
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    /**
     * @brief calculates the value of a specified variable
     * @param rValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(Parameters& rValues,
                           const Variable<double>& rThisVariable,
                           double& rValue) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Small Strain Isotropic Damage 3D constitutive law\n";
    };

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}

    ///@name Protected member Variables
    ///@{
    bool mInelasticFlag; /// Flags when in inelastic regime
    double mStrainVariable;
    ///@}

    ///@name Protected Operators
    ///@{
    ///@}

    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the stress and constitutive tensor
     * @param rValues The norm of the deviation stress
     * @param rStrainVariable
     */
    virtual void CalculateStressResponse(
            ConstitutiveLaw::Parameters &rValues,
            double &rStrainVariable);

    /**
     * @brief This method computes the positive stress vector, which in the traction-only model, is different from the stress vector.
     * @param rStressVectorPos
     * @param rStressVector
     */
    virtual void ComputePositiveStressVector(
            Vector& rStressVectorPos,
            Vector& rStressVector);

    /**
     * @brief Computes H(r), the hardening module value as a function of the strain variable
     * @param StrainVariable The properties of the material
     * @param rMaterialProperties The elastic tensor/matrix to be computed
     */
    double EvaluateHardeningModulus(
            double StrainVariable,
            const Properties &rMaterialProperties);

    /**
     * @brief Computes q(r), the hardening law value as a function of the strain variable
     * @param StrainVariable The properties of the material
     * @param rMaterialProperties The elastic tensor/matrix to be computed
     */
    double EvaluateHardeningLaw(
            double StrainVariable,
            const Properties &rMaterialProperties);

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}

    ///@name Member Variables
    ///@{

    ///@}

    ///@name Private Operators
    ///@{

    ///@}

    ///@name Private Operations
    ///@{

    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // class LinearIsotropicDamage3DLaw
} // namespace Kratos
#endif
