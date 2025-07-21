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
//  Main author:    Alejandro Cornejo/Sergio Jiménez
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "generic_small_strain_high_cycle_fatigue_law.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}GenericSmallStrainHighCycleFatiguePlaneStressLaw
///@name Kratos Classes
///@{
/**
 * @class GenericSmallStrainHighCycleFatiguePlaneStressLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the base class which define all the constitutive laws for damage in plane stress small deformation
 * @details This class considers a constitutive law integrator as an intermediate utility to compute the damage
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author Alejandro Cornejo
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericSmallStrainHighCycleFatiguePlaneStressLaw
    : public  GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base class
    typedef GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{
    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainHighCycleFatiguePlaneStressLaw);
    /**
    * Default constructor.
    */
    GenericSmallStrainHighCycleFatiguePlaneStressLaw()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainHighCycleFatiguePlaneStressLaw<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainHighCycleFatiguePlaneStressLaw(const GenericSmallStrainHighCycleFatiguePlaneStressLaw &rOther)
        : BaseType(rOther)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainHighCycleFatiguePlaneStressLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Returns the value of a specified variable (matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters &rParameterValues,
        const Variable<Matrix> &rThisVariable,
        Matrix &rValue) override;

    ///@}
    ///@name Operations
    ///@{

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    ///@}

};

} // namespace Kratos
