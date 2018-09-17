// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi

#if !defined(KRATOS_LINEAR_ISOTROPIC_DAMAGE_PLANE_STRAIN_2D_H_INCLUDED)
#define KRATOS_LINEAR_ISOTROPIC_DAMAGE_PLANE_STRAIN_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "linear_isotropic_damage_3D_law.h"

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
 * @class LinearIsotropicDamagePlaneStrain2D
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a damage with hardening/softening constitutive law in 2D (Plane Strain)
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * - INFINITY_YIELD_STRESS
 * - ISOTROPIC_HARDENING_MODULUS
 * @warning Valid for small strains
 * @note
 * @author Marcelo Raschi
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearIsotropicDamagePlaneStrain2D
    : public LinearIsotropicDamage3D

{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo ProcessInfoType;
    typedef LinearIsotropicDamage3D BaseType;
    typedef std::size_t SizeType;

    // Counted pointer of LinearIsotropicDamagePlaneStrain2DLaw
    KRATOS_CLASS_POINTER_DEFINITION(LinearIsotropicDamagePlaneStrain2D);

    /**
     * @brief Default constructor.
     */
    LinearIsotropicDamagePlaneStrain2D();

    /**
     * @brief Copy constructor.
     */
    LinearIsotropicDamagePlaneStrain2D(const LinearIsotropicDamagePlaneStrain2D& rOther);

    /**
     * @brief Destructor.
     */
    ~LinearIsotropicDamagePlaneStrain2D() override;

    /**
    * @brief Clone function
    * @return a pointer to a new instance of this constitutive law
    */
    ConstitutiveLaw::Pointer Clone() const override;

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
        return 3;
    };

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Linear Isotropic Damage Plane Strain 2D constitutive law\n";
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
    void CalculateConstitutiveMatrix(
            Matrix &constitutiveMatrix,
            const Properties &rMaterialProperties
            ) override;
    ///@}

    ///@name Protected Operations
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

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // class LinearIsotropicDamagePlaneStrain2DLaw
} // namespace Kratos
#endif
