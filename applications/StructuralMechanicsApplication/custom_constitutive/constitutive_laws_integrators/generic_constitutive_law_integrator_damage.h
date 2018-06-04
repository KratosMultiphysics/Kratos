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
#define  KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_DAMAGE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"

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
 * @brief: This object integrates the predictive stress using the damage theory by means of 
 * linear/exponential softening
 * @details
 * @tparam TYieldSurfaceType
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TYieldSurfaceType, class TVoigtSize>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegratorDamage
{
public:
    ///@name Type Definitions
    ///@{


    /// The type of yield surface
    typedef typename TYieldSurfaceType YieldSurfaceType;

    /// The type of plastic potential 
    typedef typename YieldSurfaceType::TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegratorDamage
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegratorDamage);

    /// Initialization constructor
    GenericConstitutiveLawIntegratorDamage()
    {
    }

    /// Copy constructor
    GenericConstitutiveLawIntegratorDamage(GenericConstitutiveLawIntegratorDamage const& rOther)
    {
    }

    /// Assignment operator
    GenericConstitutiveLawIntegratorDamage& operator=(GenericConstitutiveLawIntegratorDamage const& rOther)
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

    static void IntegrateStressVector(
        Vector& PredictiveStressVector,
        const double UniaxialStress,
        double& Damage,
        double& Threshold,
        const Properties& rMaterialProperties,
        const double CharacteristicLength
    )
    {
        const std::string& SofteningType = rMaterialProperties[SOFTENING_TYPE];
        double DamageParameter;
        TYieldSurfaceType::CalculateDamageParameter(rMaterialProperties, DamageParameter, CharacteristicLength);

        switch(SofteningType)
        {
            case "Linear":
                CalculateLinearDamage(UniaxialStress, Threshold, DamageParameter,
                    CharacteristicLength, rMaterialProperties, Damage);
                break;

            case "Exponential":
                CalculateExponentialDamage(UniaxialStress, Threshold, DamageParameter,
                    CharacteristicLength, rMaterialProperties, Damage);
                break;
            

            // Add more...
        }

        PredictiveStressVector *= (1.0 - Damage);
        Threshold = UniaxialStress;
    }

    static void CalculateExponentialDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        const Properties& rMaterialProperties,
        double& Damage
    )
    {
        double InitialThreshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, InitialThreshold);

        Damage = 1 - (InitialThreshold / UniaxialStress)*std::exp(DamageParameter*(1 - UniaxialStress / InitialThreshold)); 
    }

    static void CalculateLinearDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        const Properties& rMaterialProperties,
        double& Damage
    )
    {
        // TODO
    }


    static void CalculateDeviatorVector(const Vector& StressVector, Vector& rDeviator, const double rJ2)
    {
        rDeviator = StressVector;
        const double I1 = StressVector[0] + StressVector[1] + StressVector[2];
        const double Pmean = I1 / 3.0;

        rDeviator[0] -= Pmean;
        rDeviator[1] -= Pmean;
        rDeviator[2] -= Pmean;

        rJ2 = 0.5*(rDeviator[0]*rDeviator[0] + Deviator[1]*rDeviator[1] + rDeviator[2]*rDeviator[2]) +
            (rDeviator[3]*rDeviator[3] + rDeviator[4]*rDeviator[4] + rDeviator[5]*rDeviator[5]);
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

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
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

}// namespace Kratos.
#endif
