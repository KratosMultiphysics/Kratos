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
//  Main author:     Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "generic_finite_strain_kinematic_plasticity.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}GenericFiniteStrainKinematicPlasticityPlaneStress
///@name Kratos Classes
///@{
/**
 * @class GenericFiniteStrainKinematicPlasticityPlaneStress
 * @ingroup StructuralMechanicsApplication
 * @brief This class is a kinematic plasticity derived class in which we implement the plane stress condition
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author Alejandro Cornejo
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericFiniteStrainKinematicPlasticityPlaneStress
    : public  GenericFiniteStrainKinematicPlasticity<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base class
    using BaseType = GenericFiniteStrainKinematicPlasticity<TConstLawIntegratorType>;

    ///@}
    ///@name Life Cycle
    ///@{
    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericFiniteStrainKinematicPlasticityPlaneStress);

    /**
    * Default constructor.
    */
    GenericFiniteStrainKinematicPlasticityPlaneStress()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericFiniteStrainKinematicPlasticityPlaneStress<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericFiniteStrainKinematicPlasticityPlaneStress(const GenericFiniteStrainKinematicPlasticityPlaneStress &rOther)
        : BaseType(rOther)
    {
    }

    /**
    * Destructor.
    */
    ~GenericFiniteStrainKinematicPlasticityPlaneStress() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
    * @brief It calculates the constitutive matrix rConstitutiveMatrix
    * @param rConstitutiveMatrix The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues
        ) override
    {
        const auto &r_props = rValues.GetMaterialProperties();
        const double E = r_props[YOUNG_MODULUS];
        const double NU = r_props[POISSON_RATIO];
        ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStress(rConstitutiveMatrix, E, NU);
    }

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
