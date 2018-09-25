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
#include "includes/checks.h"
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"

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
 * @class ModifiedMohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Modified Mohr-Coulumb theory
 * @details The Mohr–Coulomb yield (failure) criterion is similar to the Tresca criterion, with additional provisions for materials with different tensile and compressive yield strengths. This model is often used to model concrete, soil or granular materials. This is a modified version of the criteria
 * The yield surface requires the definition of the following properties:
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @see https://en.wikipedia.org/wiki/Mohr%E2%80%93Coulomb_theory
 * @tparam TPlasticPotentialType The plastic potential considered
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

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;
    
    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;
    
    /// Counted pointer of ModifiedMohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ModifiedMohrCoulombYieldSurface);

    /// The machine precision zero tolerance
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
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rEquivalentStress The effective stress or equivalent uniaxial stress is a scalar. It is an invariant value which measures the “intensity” of a 3D stress state.
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStress(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double yield_compression = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("ModifiedMohrCoulombYieldSurface") << "Friction Angle not defined, assumed equal to 32 deg " << std::endl;
        }

        double theta;
        const double R = std::abs(yield_compression / yield_tension);
        const double Rmorh = std::pow(std::tan((Globals::Pi / 4.0) + friction_angle / 2.0), 2);
        const double alpha_r = R / Rmorh;
        const double sin_phi = std::sin(friction_angle);

        double I1, J2, J3;
        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        array_1d<double, VoigtSize> deviator(6, 0.0);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(deviator, J3);

        const double K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
        const double K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
        const double K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);

        // Check Modified Mohr-Coulomb criterion
        if (I1 == 0.0) {
            rEquivalentStress = 0.0;
        } else {
            ConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, theta);
            rEquivalentStress = (2.0 * std::tan(Globals::Pi * 0.25 + friction_angle * 0.5) / std::cos(friction_angle)) * ((I1 * K3 / 3.0) +
                        std::sqrt(J2) * (K1 * std::cos(theta) - K2 * std::sin(theta) * sin_phi / std::sqrt(3.0)));
        }
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThreshold(ConstitutiveLaw::Parameters& rValues, double& rThreshold)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const double yield_compression = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        rThreshold = std::abs(yield_compression);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviator
     * @param rg The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& GFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, GFlux, rValues);
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
        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double yield_compression = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        const double n = yield_compression / yield_tension;

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (fracture_energy * n * n * young_modulus / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
            KRATOS_DEBUG_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else { // linear
            rAParameter = -std::pow(yield_compression, 2) / (2.0 * young_modulus * fracture_energy * n * n / CharacteristicLength);
        }
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rPredictiveStressVector The stress vector
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviator
     * @param rFFlux The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateYieldSurfaceDerivative(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rFFlux,
        ConstitutiveLaw::Parameters& rValues)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;

        ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        ConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        const double checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        double c1, c2, c3;
        const double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0;
        const double sin_phi = std::sin(friction_angle);
        const double cons_phi = std::cos(friction_angle);
        const double sin_theta = std::sin(lode_angle);
        const double cos_theta = std::cos(lode_angle);
        const double cos_3theta = std::cos(3.0 * lode_angle);
        const double tan_theta = std::tan(lode_angle);
        const double tan_3theta = std::tan(3.0 * lode_angle);
        const double Root3 = std::sqrt(3.0);

        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double compr_yield = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double tens_yield = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        const double n = compr_yield / tens_yield;

        const double dilatancy = r_material_properties[DILATANCY_ANGLE] * Globals::Pi / 180.0;
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

        if (checker < 29.0) {
            c2 = cos_theta * CFL * (K1 * (1 + tan_theta * tan_3theta) + K2 * sin_phi * (tan_3theta - tan_theta) / Root3);
            c3 = CFL * (K1 * Root3 * sin_theta + K2 * sin_phi * cos_theta) / (2.0 * J2 * cos_3theta);
        } else {
            c3 = 0.0;
            double aux = 1.0;
            if (lode_angle > tolerance)
                aux = -1.0;
            c2 = 0.5 * CFL * (K1 * Root3 + aux * K2 * sin_phi / Root3);
        }
        noalias(rFFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
    }

    /**
     * @brief This method defines the check to be performed in the yield surface
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_TENSION);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_COMPRESSION);
        KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);

        if (!rMaterialProperties.Has(YIELD_STRESS)) {
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS_TENSION)) << "YIELD_STRESS_TENSION is not a defined value" << std::endl;
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION)) << "YIELD_STRESS_COMPRESSION is not a defined value" << std::endl;
            
            const double yield_compression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
            const double yield_tension = rMaterialProperties[YIELD_STRESS_TENSION];

            KRATOS_ERROR_IF(yield_compression < tolerance) << "Yield stress in compression almost zero or negative, include YIELD_STRESS_COMPRESSION in definition";
            KRATOS_ERROR_IF(yield_tension < tolerance) << "Yield stress in tension almost zero or negative, include YIELD_STRESS_TENSION in definition";
        } else {
            const double yield_stress = rMaterialProperties[YIELD_STRESS];

            KRATOS_ERROR_IF(yield_stress < tolerance) << "Yield stress almost zero or negative, include YIELD_STRESS in definition";
        }
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRACTURE_ENERGY)) << "FRACTURE_ENERGY is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not a defined value" << std::endl;

        return TPlasticPotentialType::Check(rMaterialProperties);
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
