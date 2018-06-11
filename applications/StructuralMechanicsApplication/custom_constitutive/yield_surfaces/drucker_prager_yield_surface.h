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

#if !defined(KRATOS_DRUCKER_PRAGER_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_DRUCKER_PRAGER_YIELD_SURFACE_H_INCLUDED

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
 * @class DruckerPragerYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
//template <class TPlasticPotentialType , std::size_t TVoigtSize>
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DruckerPragerYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of DruckerPragerYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(DruckerPragerYieldSurface);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    DruckerPragerYieldSurface()
    {
    }

    /// Copy constructor
    DruckerPragerYieldSurface(DruckerPragerYieldSurface const& rOther)
    {
    }

    /// Assignment operator
    DruckerPragerYieldSurface& operator=(DruckerPragerYieldSurface const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~DruckerPragerYieldSurface() {};

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
        double friction_angle = rMaterialProperties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
        const double SinPhi = std::sin(friction_angle);
        const double Root3 = std::sqrt(3.0);

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("DruckerPragerYieldSurface") << "Friction Angle not defined, assumed equal to 32 " << std::endl;
        }

        double I1, J2;
        ConstitutiveLawUtilities::CalculateI1Invariant(StressVector, I1);
        Vector Deviator = ZeroVector(TVoigtSize);
        ConstitutiveLawUtilities::CalculateJ2Invariant(StressVector,I1, Deviator, J2);

        if (I1 == 0.0) { rEqStress = 0; }
        else {
            const double CFL = -Root3*(3.0 - SinPhi) / (3.0 * SinPhi - 3.0);
            const double TEN0 = 6.0 * I1*SinPhi / (Root3*(3.0 - SinPhi)) + std::sqrt(J2);
            rEqStress = std::abs(CFL*TEN0);
        }
    }

    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        const double YieldTension = rMaterialProperties[YIELD_STRESS_TENSION];
        const double friction_angle = rMaterialProperties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
        const double SinPhi = std::sin(friction_angle);
        
        rThreshold = std::abs(YieldTension*(3.0 + SinPhi) / (3.0*SinPhi - 3.0));
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
        c3 = 0.0;

        const double FrictionAngle = rMaterialProperties[INTERNAL_FRICTION_ANGLE];
        const double SinPhi    = std::sin(FrictionAngle);
        const double Root3     = std::sqrt(3.0);

        const double CFL = -Root3*(3.0-SinPhi) / (3.0*SinPhi-3.0);
        c1 = CFL*2.0*SinPhi / (Root3*(3.0-SinPhi));
        c2 = CFL;

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

}; // Class DruckerPragerYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
