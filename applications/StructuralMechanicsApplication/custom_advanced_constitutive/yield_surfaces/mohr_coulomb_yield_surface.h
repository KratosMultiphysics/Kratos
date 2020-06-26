// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED

// System includes

// Project includes
#include "includes/checks.h"
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"

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

///@}
///@name Kratos Classes
///@{
/**
 * @class MohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Von-Mises theory
 * @details The Mohr–Coulomb failure surface is a cone with a hexagonal cross section in deviatoric stress space
 * The yield surface requires the definition of the following properties:
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * - COHESION: Is the intercept of the failure envelope with the tau axis
 * @see https://en.wikipedia.org/wiki/Mohr%E2%80%93Coulomb_theory
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class MohrCoulombYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// Counted pointer of MohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    MohrCoulombYieldSurface()
    {
    }

    /// Copy constructor
    MohrCoulombYieldSurface(MohrCoulombYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    MohrCoulombYieldSurface &operator=(MohrCoulombYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~MohrCoulombYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStress(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double I1, J2, J3, lode_angle;
        array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);

        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(deviator, J3);
        ConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0;

        rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(friction_angle) / std::sqrt(3.0)) * std::sqrt(J2) +
            I1 * std::sin(friction_angle) / 3.0;
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThreshold(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double cohesion = r_material_properties[COHESION];
        const double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0;

        rThreshold = cohesion * std::cos(friction_angle);
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameter(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double fracture_energy = r_material_properties[FRACTURE_ENERGY];
        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        double equivalent_yield;
        GetInitialUniaxialThreshold(rValues, equivalent_yield);
        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(equivalent_yield, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture Energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else { // linear
            rAParameter = -std::pow(equivalent_yield, 2) / (2.0 * young_modulus * fracture_energy / CharacteristicLength);
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rDerivativePlasticPotential,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rFFlux The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateYieldSurfaceDerivative(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rFFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;
		const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0;

        ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        ConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        double c1, c3, c2;
		double checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        if (std::abs(checker) < 29.0) { // If it is not the edge
            c1 = std::sin(friction_angle);
            c3 = (std::sqrt(3.0) * std::sin(lode_angle) + std::sin(friction_angle) * std::cos(lode_angle)) /
                (2.0 * J2 * std::cos(3.0 * lode_angle));
            c2 = 0.5 * std::cos(lode_angle)*(1.0 + std::tan(lode_angle) * std::sin(3.0 * lode_angle) +
                std::sin(friction_angle) * (std::tan(3.0 * lode_angle) - std::tan(lode_angle)) / std::sqrt(3.0));
        } else { // smoothing with drucker-praguer
            c1 = 3.0 * (2.0 * std::sin(friction_angle) / (std::sqrt(3.0) * (3.0 - std::sin(friction_angle))));
            c2 = 1.0;
            c3 = 0.0;
        }

        noalias(rFFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
    }

    /**
     * @brief This method defines the check to be performed in the yield surface
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(COHESION);
        KRATOS_CHECK_VARIABLE_KEY(FRICTION_ANGLE);
        KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(COHESION)) << "COHESION is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRICTION_ANGLE)) << "FRICTION_ANGLE is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRACTURE_ENERGY)) << "FRACTURE_ENERGY is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS)) << "YIELD_STRESS is not a defined value" << std::endl;

        return TPlasticPotentialType::Check(rMaterialProperties);
    }

    /**
     * @brief This method returns true if the yield surfacecompares with the tension tield stress
     */
    static bool IsWorkingWithTensionThreshold()
    {
        return true;
    }

    /**
     * @brief This method returns the scaling factor of the yield surface surfacecompares with the tension tield stress
     */
    static double GetScaleFactorTension(const Properties& rMaterialProperties)
    {
        return 1.0;
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

    ///@}

}; // Class MohrCoulombYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
