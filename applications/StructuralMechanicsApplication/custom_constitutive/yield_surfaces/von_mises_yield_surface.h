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

#if !defined(KRATOS_VON_MISES_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_VON_MISES_YIELD_SURFACE_H_INCLUDED

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
 * @class VonMisesYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType 
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VonMisesYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of VonMisesYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( VonMisesYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    VonMisesYieldSurface()
    {
    }

    /// Copy constructor
    VonMisesYieldSurface(VonMisesYieldSurface const& rOther)
    {
    }

    /// Assignment operator
    VonMisesYieldSurface& operator=(VonMisesYieldSurface const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~VonMisesYieldSurface() {};

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
        double I1, J2;
        Vector Deviator = ZeroVector(TVoigtSize);

        ConstitutiveLawUtilities::CalculateI1Invariant(StressVector, I1);
        ConstitutiveLawUtilities::CalculateJ2Invariant(StressVector, I1, Deviator, J2);

        rEqStress = std::sqrt(3.0*J2);
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

    /*
    This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
    As:            DF/DS = c1*V1 + c2*V2 + c3*V3
    */
    static void CalculateYieldSurfaceDerivative(
        const Vector& StressVector, 
        const Vector& Deviator,
        const double J2, 
        Vector& rFFlux,
        const Properties& rMaterialProperties
        )
    {
        Vector FirstVector, SecondVector, ThirdVector;

        ConstitutiveLawUtilities::CalculateFirstVector(FirstVector);
        ConstitutiveLawUtilities::CalculateSecondVector(Deviator, J2, SecondVector);
        ConstitutiveLawUtilities::CalculateThirdVector(Deviator, J2, ThirdVector);

        double c1, c2, c3;
        c1 = 0.0;
        c2 = std::sqrt(3.0);
        c3 = 0.0;

        noalias(rFFlux) = c1*FirstVector + c2*SecondVector + c3*ThirdVector;
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

}; // Class VonMisesYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
