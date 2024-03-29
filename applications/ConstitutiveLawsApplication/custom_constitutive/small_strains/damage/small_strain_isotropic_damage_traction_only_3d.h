// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//  Collaborator:

#pragma once

// System includes

// External includes

// Project includes
#include "small_strain_isotropic_damage_3d.h"

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
 * @class LinearIsotropicDamageTractionOnly3D
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a damage with hardening/softening constitutive law in 3D
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
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SmallStrainIsotropicDamageTractionOnly3D
    : public SmallStrainIsotropicDamage3D

{
public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    // Counted pointer of LinearIsotropicDamageTractionOnly3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainIsotropicDamageTractionOnly3D);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainIsotropicDamageTractionOnly3D();

    /**
     * @brief Default constructor.
     */
    SmallStrainIsotropicDamageTractionOnly3D(const SmallStrainIsotropicDamageTractionOnly3D& rOther);

    /**
     * @brief Default constructor.
     */
    ~SmallStrainIsotropicDamageTractionOnly3D() override;

    /**
     * @brief Clone function
     * @return A pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

//    /**
//     * @brief This method computes the stress and constitutive tensor
//     * @param rValues The norm of the deviation stress
//     * @param rStrainVariable
//     */
//    void CalculateStressResponse(
//            ConstitutiveLaw::Parameters& rValues,
//            double& rStrainVariable) override;

    /**
     * @brief This method computes the positive stress vector, which in the traction-only model, is different from the stress vector.
     * @param rStressVectorPos
     * @param rStressVector
     */
    void ComputePositiveStressVector(
            Vector& rStressVectorPos,
            Vector& rStressVector) override;

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
