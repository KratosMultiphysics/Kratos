// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//

#if !defined(KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_DAMAGE_H_INCLUDED)
#define KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_DAMAGE_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"

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
 * @class GenericConstitutiveLawIntegratorDamage
 * @ingroup StructuralMechanicsApplication
 * @brief: This object integrates the predictive stress using the isotropic damage theory by means of 
 * linear/exponential softening.
 * @details
 * @tparam TYieldSurfaceType
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TYieldSurfaceType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegratorDamage
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of yield surface
    typedef TYieldSurfaceType YieldSurfaceType;

    /// The type of plastic potential
    typedef typename YieldSurfaceType::PlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegratorDamage
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegratorDamage);

    /// Initialization constructor
    GenericConstitutiveLawIntegratorDamage()
    {
    }

    /// Copy constructor
    GenericConstitutiveLawIntegratorDamage(GenericConstitutiveLawIntegratorDamage const &rOther)
    {
    }

    /// Assignment operator
    GenericConstitutiveLawIntegratorDamage &operator=(GenericConstitutiveLawIntegratorDamage const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericConstitutiveLawIntegratorDamage()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method integrates the predictive stress vector with the CL using linear or exponential softening
     * @param PredictiveStressVector The predictive stress vector
     * @param UniaxialStress The equivalent uniaxial stress 
     * @param Damage The internal variable of the damage model
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void IntegrateStressVector(
        Vector &rPredictiveStressVector,
        const double UniaxialStress,
        double &rDamage,
        double &rThreshold,
        const Properties &rMaterialProperties,
        const double CharacteristicLength)
    {
        const int softening_type = rMaterialProperties[SOFTENING_TYPE];
        double DamageParameter;
        TYieldSurfaceType::CalculateDamageParameter(rMaterialProperties, DamageParameter, CharacteristicLength);

        switch (softening_type)
        {
        case static_cast<int>(SofteningType::Linear):
            CalculateLinearDamage(UniaxialStress, rThreshold, DamageParameter,
                                  CharacteristicLength, rMaterialProperties, rDamage);
            break;
        case static_cast<int>(SofteningType::Exponential):
            CalculateExponentialDamage(UniaxialStress, rThreshold, DamageParameter,
                                       CharacteristicLength, rMaterialProperties, rDamage);
            break;
        default:
            KRATOS_ERROR << "SOFTENING_TYPE not defined or wrong..." << softening_type << std::endl;
            break;
        }
        rPredictiveStressVector *= (1.0 - rDamage);
    }

    /**
     * @brief This computes the damage variable according to exponential softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param Damage The internal variable of the damage model
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateExponentialDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        const Properties &rMaterialProperties,
        double &Damage)
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);
        Damage = 1.0 - (initial_threshold / UniaxialStress) * std::exp(DamageParameter * (1.0 - UniaxialStress / initial_threshold));
    }

    /**
     * @brief This computes the damage variable according to linear softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param Damage The internal variable of the damage model
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateLinearDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        const Properties &rMaterialProperties,
        double &Damage)
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);
        Damage = (1.0 - initial_threshold / UniaxialStress) / (1.0 + DamageParameter);
    }

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

    void save(Serializer &rSerializer) const
    {
    }

    void load(Serializer &rSerializer)
    {
    }

    ///@}

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
