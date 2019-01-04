// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Fernando Rastellini
//  Collaborators:   Riccardo Rossi
//                   Vicente Mataix Ferrandiz

#if !defined(KRATOS_PLASTICITY_ISOTROPIC_KINEMATIC_J2_H_INCLUDED)
#define KRATOS_PLASTICITY_ISOTROPIC_KINEMATIC_J2_H_INCLUDED

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
 * @class PlasticityIsotropicKinematicJ2
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a Simo J2 plasticity CL with Isotropic & Kinematic Hardening in 3D
 * @details This constitutive law is defined by the following parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * - ISOTROPIC_HARDENING_MODULUS
 * @warning TODO: Kinematic hardening
 * @warning Only small strains (for the moment...)
 * @author Fernando Rastellini
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PlasticityIsotropicKinematicJ2
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /// The Voigt size definition
    static constexpr SizeType VoigtSize = 6;   //3D only (for the moment)

    /// The dimension definition
    static constexpr SizeType Dimension = VoigtSize == 6 ? 3 : 2;

    /// The strain/stress vector type definition
    typedef array_1d<double, VoigtSize> BoundedArrayType;

    /// The constitutive matrix type definition
    typedef Matrix MatrixType;

    // Counted pointer of PlasticityIsotropicKinematicJ2
    KRATOS_CLASS_POINTER_DEFINITION(PlasticityIsotropicKinematicJ2);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    PlasticityIsotropicKinematicJ2();

    /**
     * @brief Copy constructor.
     */
    PlasticityIsotropicKinematicJ2(const PlasticityIsotropicKinematicJ2& rOther);

    /**
     * @brief Destructor.
     */
    ~PlasticityIsotropicKinematicJ2() override;

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
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return VoigtSize;
    };

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     */
    double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     */
    void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @warning This function is deprecated.
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
     * @warning This function is deprecated.
     */
    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues,
                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Initializes the material response in terms of 1st Piola-Kirchhoff stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initializes the material response in terms of 2nd Piola-Kirchhoff stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initializes the material response in terms of Kirchhoff stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initializes the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

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

    BoundedArrayType mPlasticStrain;           /// The last converged plastic strain
    double           mEquivalentPlasticStrain; /// The last converged equivalent plastic strain

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the yield function value
     * @param NormDeviatoricStress The euclidean norm of the deviatoric stress
     * @param rMaterialProperties The material properties considered
     * @return The yield function value (stress units)
     */
    double YieldFunction(
        const double      NormDeviatoricStress,
        const Properties& rMaterialProperties
        );

    /**
     * @brief This method computes the stress and constitutive tensor for VoigtSize = 6
     * @param rValues The norm of the deviation stress
     * @param rPlasticStrain
     * @param rAccumulatedPlasticStrain
     */
    void CalculateResponse6(
        ConstitutiveLaw::Parameters& rValues,
        BoundedArrayType&            rPlasticStrain,
        double&                      rAccumulatedPlasticStrain
        );

    /**
     * @brief This method computes the constitutive tangent matrix for VoigtSize = 6
     * @param DeltaGamma The increment on the Gamma parameter
     * @param NormStressTrial The norm of the stress trial
     * @param YieldFunctionNormalVector The yield function normal vector
     * @param rMaterialProperties The properties of the material
     * @param rTangentTensor The tangent tensor/matrix to be computed
     */
    void CalculateTangentTensor6(
        const double            DeltaGamma,
        const double            NormStressTrial,
        const BoundedArrayType& rYieldFunctionNormalVector,
        const Properties&       rMaterialProperties,
        MatrixType&             rTangentTensor
        );

    /**
     * @brief This method computes the elastic tensor for VoigtSize = 6
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     * @param rMaterialProperties The properties of the material
     */
    void CalculateElasticMatrix6(
        const Properties& rMaterialProperties,
        MatrixType&       rElasticityTensor
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

}; // Class PlasticityIsotropicKinematicJ2
} // namespace Kratos.
#endif
