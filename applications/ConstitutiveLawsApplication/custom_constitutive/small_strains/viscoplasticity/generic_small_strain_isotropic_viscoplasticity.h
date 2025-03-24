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
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"

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
 * @class GenericSmallStrainIsotropicViscoPlasticity
 * @ingroup ConstitutiveLawsApp
 * @brief 
 * @author Alejandro Cornejo
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericSmallStrainIsotropicViscoPlasticity
    : public GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    // The size type definition
    using SizeType = std::size_t;

    using ConstLawIntegratorType = TConstLawIntegratorType;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TConstLawIntegratorType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TConstLawIntegratorType::VoigtSize;

    /// Definition of the base class
    using BaseType = GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>;

    /// Counted pointer of GenericSmallStrainIsotropicPlasticity
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainIsotropicViscoPlasticity);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
   GenericSmallStrainIsotropicViscoPlasticity()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
   GenericSmallStrainIsotropicViscoPlasticity(const GenericSmallStrainIsotropicViscoPlasticity &rOther)
        : BaseType(rOther)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainIsotropicViscoPlasticity() override
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
