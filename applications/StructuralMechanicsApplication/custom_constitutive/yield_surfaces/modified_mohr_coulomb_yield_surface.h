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

#if !defined(KRATOS_MODIFIED_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_MODIFIED_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED

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
 * @class ModifiedMohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ModifiedMohrCoulombYieldSurface
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of ModifiedMohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ModifiedMohrCoulombYieldSurface);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ModifiedMohrCoulombYieldSurface()
    {
    }

    /// Copy constructor
    ModifiedMohrCoulombYieldSurface(ModifiedMohrCoulombYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    ModifiedMohrCoulombYieldSurface &operator=(ModifiedMohrCoulombYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ModifiedMohrCoulombYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param StressVector The stress vector
     * @param StrainVector The StrainVector vector
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStress(
        const Vector &StressVector,
        const Vector &StrainVector,
        double &rEqStress,
        const Properties &rMaterialProperties)
    {
        const double yield_compression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = rMaterialProperties[YIELD_STRESS_TENSION];
        double friction_angle = rMaterialProperties[FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!

        // Check input variables
        if (friction_angle < tolerance)
        {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("ModifiedMohrCoulombYieldSurface") << "Friction Angle not defined, assumed equal to 32 deg " << std::endl;
        }

        KRATOS_ERROR_IF(yield_compression < tolerance) << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_COMPRESSION in .mdpa ";
        KRATOS_ERROR_IF(yield_tension < tolerance) << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_TENSION in .mdpa ";

        double K1, K2, K3, Rmorh, R, alpha_r, theta;
        R = std::abs(yield_compression / yield_tension);
        Rmorh = std::pow(std::tan((Globals::Pi / 4.0) + friction_angle / 2.0), 2);
        alpha_r = R / Rmorh;
        double sin_phi = std::sin(friction_angle);

        double I1, J2, J3;
        ConstitutiveLawUtilities::CalculateI1Invariant(StressVector, I1);
        Vector Deviator = ZeroVector(6);
        ConstitutiveLawUtilities::CalculateJ2Invariant(StressVector, I1, Deviator, J2);
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);

        K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
        K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
        K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);

        // Check Modified Mohr-Coulomb criterion
        if (I1 == 0.0)
        {
            rEqStress = 0.0;
        }
        else
        {
            ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, theta);
            rEqStress = (2.0 * std::tan(Globals::Pi * 0.25 + friction_angle * 0.5) / std::cos(friction_angle)) * ((I1 * K3 / 3.0) +
                        std::sqrt(J2) * (K1 * std::cos(theta) - K2 * std::sin(theta) * sin_phi / std::sqrt(3.0)));
        }
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rMaterialProperties The material properties
     */
    static void GetInitialUniaxialThreshold(const Properties &rMaterialProperties, double &rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_COMPRESSION]);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rg The derivative of the plastic potential
     * @param rMaterialProperties The material properties
     */
    static void CalculatePlasticPotentialDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &GFlux,
        const Properties &rMaterialProperties)
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, GFlux, rMaterialProperties);
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param AParameter The damage parameter
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameter(
        const Properties &rMaterialProperties,
        double &AParameter,
        const double CharacteristicLength)
    {
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];
        const double E = rMaterialProperties[YOUNG_MODULUS];
        const double sigma_c = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double sigma_t = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = sigma_c / sigma_t;

        if (rMaterialProperties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential))
        {
            AParameter = 1.00 / (Gf * n * n * E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
            KRATOS_ERROR_IF(AParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        }
        else
        { // linear
            AParameter = -std::pow(sigma_c, 2) / (2.0 * E * Gf * n * n / CharacteristicLength);
        }
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rFFlux The derivative of the yield surface
     * @param rMaterialProperties The material properties
     */
    static void CalculateYieldSurfaceDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &rFFlux,
        const Properties &rMaterialProperties)
    {
        Vector FirstVector, SecondVector, ThirdVector;

        ConstitutiveLawUtilities::CalculateFirstVector(FirstVector);
        ConstitutiveLawUtilities::CalculateSecondVector(Deviator, J2, SecondVector);
        ConstitutiveLawUtilities::CalculateThirdVector(Deviator, J2, ThirdVector);

        double J3, lode_angle;
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, lode_angle);

        const double Checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        double c1, c2, c3;
        const double friction_angle = rMaterialProperties[FRICTION_ANGLE] * Globals::Pi / 180.0;
        const double sin_phi = std::sin(friction_angle);
        const double cons_phi = std::cos(friction_angle);
        const double sin_theta = std::sin(lode_angle);
        const double cos_theta = std::cos(lode_angle);
        const double cos_3theta = std::cos(3.0 * lode_angle);
        const double tan_theta = std::tan(lode_angle);
        const double tan_3theta = std::tan(3.0 * lode_angle);
        const double Root3 = std::sqrt(3.0);

        const double compr_yield = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double tens_yield = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = compr_yield / tens_yield;

        const double dilatancy = rMaterialProperties[DILATANCY_ANGLE] * Globals::Pi / 180.0;
        ;
        const double angle_phi = (Globals::Pi * 0.25) + dilatancy * 0.5;
        const double alpha = n / (std::tan(angle_phi) * std::tan(angle_phi));

        const double CFL = 2.0 * std::tan(angle_phi) / cons_phi;

        const double K1 = 0.5 * (1.0 + alpha) - 0.5 * (1.0 - alpha) * sin_phi;
        const double K2 = 0.5 * (1.0 + alpha) - 0.5 * (1.0 - alpha) / sin_phi;
        const double K3 = 0.5 * (1.0 + alpha) * sin_phi - 0.5 * (1.0 - alpha);

        if (std::abs(sin_phi) > tolerance)
            c1 = CFL * K3 / 3.0;
        else
            c1 = 0.0; // check

        if (Checker < 29.0)
        {
            c2 = cos_theta * CFL * (K1 * (1 + tan_theta * tan_3theta) + K2 * sin_phi * (tan_3theta - tan_theta) / Root3);
            c3 = CFL * (K1 * Root3 * sin_theta + K2 * sin_phi * cos_theta) / (2.0 * J2 * cos_3theta);
        }
        else
        {
            c3 = 0.0;
            double aux = 1.0;
            if (lode_angle > tolerance)
                aux = -1.0;
            c2 = 0.5 * CFL * (K1 * Root3 + aux * K2 * sin_phi / Root3);
        }
        noalias(rFFlux) = c1 * FirstVector + c2 * SecondVector + c3 * ThirdVector;
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

}; // Class ModifiedMohrCoulombYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
