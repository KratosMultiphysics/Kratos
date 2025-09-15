// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// External includes
// Project includes
#include "custom_constitutive/thermal/small_strains/damage/generic_small_strain_thermal_isotropic_damage.h"

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
 * @class GenericSmallStrainThermalIsotropicDamage
 * @ingroup ConstitutiveLawsApp
 * @brief This class derives from the Isotropic damage CL and adds thermal effects (material properties affectation and internal variables)
 * @details This class considers a constitutive law integrator as an intermediate utility to compute the damage. 3D and plane strain
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author Alejandro Cornejo
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericSmallStrainThermalIsotropicDamagePlaneStress
    : public GenericSmallStrainThermalIsotropicDamage<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{


    /// Definition of the base class
    using BaseType = GenericSmallStrainThermalIsotropicDamage<TConstLawIntegratorType>;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = BaseType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = BaseType::VoigtSize;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainThermalIsotropicDamagePlaneStress);

    /// Advanced and basic contitutive laws utilities for the corresponding Voigt size
    using CLutils    = ConstitutiveLawUtilities<VoigtSize>;
    using AdvCLutils = AdvancedConstitutiveLawUtilities<VoigtSize>;

    /// Bounded vector for stresses/strains
    using BoundedArrayType = array_1d<double, VoigtSize>;

    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();
    static constexpr double threshold_tolerance = 1.0e-5;


    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainThermalIsotropicDamagePlaneStress()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainThermalIsotropicDamagePlaneStress<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainThermalIsotropicDamagePlaneStress(const GenericSmallStrainThermalIsotropicDamagePlaneStress &rOther)
        : BaseType(rOther)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainThermalIsotropicDamagePlaneStress() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;


    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;


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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GenericSmallStrainThermalIsotropicDamage<TConstLawIntegratorType>)
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GenericSmallStrainThermalIsotropicDamage<TConstLawIntegratorType>)
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
