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

#if !defined(KRATOS_LINEAR_J2_PLASTIC_PLANE_STRAIN_2D_H_INCLUDED)
#define KRATOS_LINEAR_J2_PLASTIC_PLANE_STRAIN_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "linear_j2_plasticity_3d.h"

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
 * @class LinearJ2PlasticityPlaneStrain2D
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a Simo J2 plasticity constitutive law in 2D (Plane Strain)
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * - REFERENCE_HARDENING_MODULUS (kinematic hardening modulus)
 * - ISOTROPIC_HARDENING_MODULUS
 * - INFINITY_HARDENING_MODULUS (saturation hardening modulus)
 * - HARDENING_EXPONENT
 * @warning Valid for small strains, linear quadrilaterals
 * @note Strain size is 4 (xx, yy, zz, xy).
 * @note Requires B-bar element
 * @author Marcelo Raschi
 * @author Manuel Caicedo
 * @author Alfredo Huespe
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearJ2PlasticityPlaneStrain2D
    : public LinearJ2Plasticity3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw       CLBaseType;
    typedef LinearJ2Plasticity3D    BaseType;
    typedef std::size_t             SizeType;

    // Counted pointer of LinearJ2Plasticity3D
    KRATOS_CLASS_POINTER_DEFINITION(LinearJ2PlasticityPlaneStrain2D);

    /**
     * @brief Default constructor.
     */
    LinearJ2PlasticityPlaneStrain2D();

    /**
     * @brief Copy constructor.
     */
    LinearJ2PlasticityPlaneStrain2D(const LinearJ2PlasticityPlaneStrain2D& rOther);

    /**
     * @brief Destructor.
     */
    ~LinearJ2PlasticityPlaneStrain2D() override;

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
        return 2;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 4;
    };

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

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Linear J2 Plasticity Plane Strain 2D constitutive law\n";
    };

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

    /**
     * @brief This method computes the plastic potential
     * @param DeltaGamma The increment on the Gamma parameter
     * @param NormStressTrial The norm of the stress trial
     * @param YieldFunctionNormalVector The yield function normal vector
     * @param rMaterialProperties The properties of the material
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     */
    void CalculateTangentTensor(
        const double DeltaGamma,
        const double NormStressTrial,
        const Vector& YieldFunctionNormalVector,
        const Properties& rMaterialProperties,
        Matrix& rElasticityTensor
        ) override;

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     * @param rMaterialProperties The properties of the material
     */
    void CalculateElasticMatrix(
        Matrix &rElasticityTensor,
        const Properties &rMaterialProperties
        ) override;

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

}; // Class LinearJ2PlasticityPlaneStrain2D
} // namespace Kratos.
#endif
