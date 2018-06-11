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

#if !defined(KRATOS_RANKINE_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_RANKINE_YIELD_SURFACE_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"

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
 * @class RankineYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) RankineYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of RankineYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(RankineYieldSurface);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    RankineYieldSurface()
    {
    }

    /// Copy constructor
    RankineYieldSurface(RankineYieldSurface const& rOther)
    {
    }

    /// Assignment operator
    RankineYieldSurface& operator=(RankineYieldSurface const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~RankineYieldSurface() {};

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    static void CalculateEquivalentStress(  
        const Vector& StressVector,
        const Vector& StrainVector,
        double& rEqStress, 
        const Properties& rMaterialProperties
    )
    {      
        Vector PrincipalStressVector = ZeroVector(3);
        ConstitutiveLawUtilities::CalculatePrincipalStresses(PrincipalStressVector, StressVector);

        // the rEqStress is the maximum principal stress
        rEqStress = std::max(std::max(PrincipalStressVector[0],PrincipalStressVector[1]), PrincipalStressVector[2]);
    }


    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_TENSION]);
    }

    static void CalculateDamageParameter(
        const Properties& rMaterialProperties, 
        double& AParameter, 
        const double CharacteristicLength
        )
    {
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];
        const double E  = rMaterialProperties[YOUNG_MODULUS];
        const double sigma_c = rMaterialProperties[YIELD_STRESS_COMPRESSION];

        if (rMaterialProperties[SOFTENING_TYPE] == static_cast<std::size_t>(SofteningType::Exponential)) {
            AParameter = 1.00 / (Gf*E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
        } else {
            
        }
    }

    // Computes dG/dS
    static void CalculatePlasticPotentialDerivative(
        const Vector& StressVector,
        const Vector& Deviator,
        const double J2, 
        Vector& rg,
        const Properties& rMaterialProperties
    )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, rg, rMaterialProperties);
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

}; // Class RankineYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
