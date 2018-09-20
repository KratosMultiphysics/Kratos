// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Alfredo Huespe
//  Collaborator:    Vicente Mataix Ferrandiz

#if !defined(KRATOS_LINEAR_J2_PLASTIC_3D_H_INCLUDED)
#define KRATOS_LINEAR_J2_PLASTIC_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
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
 * @class LinearJ2Plasticity3D
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a Simo J2 plasticity constitutive law in 3D
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * - REFERENCE_HARDENING_MODULUS (kinematic hardening modulus)
 * - ISOTROPIC_HARDENING_MODULUS
 * - INFINITY_HARDENING_MODULUS (saturation hardening modulus)
 * - HARDENING_EXPONENT
 * @warning Valid for small strains, linear hexahedra
 * @note Requires B-bar element
 * @author Marcelo Raschi
 * @author Manuel Caicedo
 * @author Alfredo Huespe
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearJ2Plasticity3D
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    // Counted pointer of LinearJ2Plasticity3D
    KRATOS_CLASS_POINTER_DEFINITION(LinearJ2Plasticity3D);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    LinearJ2Plasticity3D();

    /**
     * @brief Copy constructor.
     */
    LinearJ2Plasticity3D(const LinearJ2Plasticity3D& rOther);

    /**
     * @brief Destructor.
     */
    ~LinearJ2Plasticity3D() override;

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
     * @brief Returns whether this constitutive Law has specified variable (bool)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

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
     * @brief To be called at the end of each solution step  (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param rCurrentProcessInfo the current ProcessInfo instance
     */
    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues,
                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
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
        rOStream << "Linear J2 Plasticity 3D constitutive law\n";
    };

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    bool mInelasticFlag; /// This flags tells if we are in a elastic or ineslastic regime
    Vector mPlasticStrain; /// The current plastic strain (one for each of the strain components)
    Vector mPlasticStrainOld; /// The previous plastic strain (one for each of the strain components)
    double mAccumulatedPlasticStrain; /// The current accumulated plastic strain
    double mAccumulatedPlasticStrainOld; /// The previous accumulated plastic strain

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the yield function
     * @param NormDeviationStress The norm of the deviation stress
     * @param rMaterialProperties The properties of the current material considered
     * @return The trial yield function (after update)
     */
    double YieldFunction(
        const double NormDeviationStress,
        const Properties& rMaterialProperties
        );

    /**
     * @brief This method computes the increment of Gamma
     * @param NormStressTrial The norm of the stress trial
     * @param rMaterialProperties The properties of the material
     * @return The increment of Gamma computed
     */
    double GetDeltaGamma(
        const double NormStressTrial,
        const Properties& rMaterialProperties
        );

    /**
     * @brief This method gets the saturation hardening parameter
     * @param rMaterialProperties The properties of the material
     * @return The saturation hardening parameter
     */
    double GetSaturationHardening(const Properties& rMaterialProperties);

    /**
     * @brief This method computes the plastic potential
     * @param rMaterialProperties The properties of the material
     * @return The plastic potential
     */
    double GetPlasticPotential(const Properties& rMaterialProperties);

    /**
     * @brief This method computes the plastic potential
     * @param DeltaGamma The increment on the Gamma parameter
     * @param NormStressTrial The norm of the stress trial
     * @param YieldFunctionNormalVector The yield function normal vector
     * @param rMaterialProperties The properties of the material
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     */
    virtual void CalculateTangentTensor(
        const double DeltaGamma,
        const double NormStressTrial,
        const Vector& YieldFunctionNormalVector,
        const Properties& rMaterialProperties,
        Matrix& rElasticityTensor
        );

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     * @param rMaterialProperties The properties of the material
     */
    virtual void CalculateElasticMatrix(
        Matrix &rElasticityTensor,
        const Properties &rMaterialProperties
        );

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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class LinearJ2Plasticity3D
} // namespace Kratos.
#endif
