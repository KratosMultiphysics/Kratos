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

#if !defined(KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED

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
 * @class TrescaYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TrescaYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of TrescaYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( TrescaYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    TrescaYieldSurface()
    {
    }

    /// Copy constructor
    TrescaYieldSurface(TrescaYieldSurface const& rOther)
    {
    }

    /// Assignment operator
    TrescaYieldSurface& operator=(TrescaYieldSurface const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~TrescaYieldSurface() {};

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
        double I1, J2, J3, LodeAngle;
        Vector Deviator = ZeroVector(TVoigtSize);

        ConstitutiveLawUtilities::CalculateI1Invariant(StressVector, I1);
        ConstitutiveLawUtilities::CalculateJ2Invariant(StressVector, I1, Deviator, J2);
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, LodeAngle);

        rEqStress = 2.0*std::cos(LodeAngle)*std::sqrt(J2);
    }

    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_TENSION]);  // TODO Check
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

        double J3, LodeAngle;
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, LodeAngle);

        const double Checker = std::abs(LodeAngle*57.29577951308);

        double c1, c2, c3;
        c1 = 0.0;

        if (Checker < 29.0) {
            c2 = 2.0*(std::cos(LodeAngle) + std::sin(LodeAngle)*std::tan(3.0*LodeAngle));
            c3 = std::sqrt(3.0)*std::sin(LodeAngle) / (J2*std::cos(3.0*LodeAngle));
        } else {
            c2 = std::sqrt(3.0);
            c3 = 0.0;
        }

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

}; // Class TrescaYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
