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

#pragma once

// System includes

// Project includes
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "thermal_von_mises_yield_surface.h"

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
 * @class ThermalMohrCoulombYieldSurface
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
 * @author Alejandro Cornejo
 */
template <class TPlasticPotentialType>
class ThermalMohrCoulombYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    using PlasticPotentialType = TPlasticPotentialType;

    using BaseType = MohrCoulombYieldSurface<TPlasticPotentialType>;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// Counted pointer of MohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ThermalMohrCoulombYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Advanced contitutive laws utilities for the corresponding Voigt size
    using AdvCLutils = AdvancedConstitutiveLawUtilities<VoigtSize>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ThermalMohrCoulombYieldSurface()
    {
    }

    /// Copy constructor
    ThermalMohrCoulombYieldSurface(ThermalMohrCoulombYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    ThermalMohrCoulombYieldSurface &operator=(ThermalMohrCoulombYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ThermalMohrCoulombYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param rStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStress(
        const array_1d<double, VoigtSize>& rStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double I1, J2, J3, lode_angle;
        array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);

        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rStressVector, I1);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rStressVector, I1, deviator, J2);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(deviator, J3);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        const double friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0;

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
        const auto& r_material_properties = rValues.GetMaterialProperties();
        double friction_angle, cohesion;
        if (rValues.IsSetShapeFunctionsValues()) { // This is needed since at Initialize level the N are not set yet...
            friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues);
            cohesion = AdvCLutils::GetMaterialPropertyThroughAccessor(COHESION, rValues);
        } else {
            const double ref_temperature = r_material_properties.Has(REFERENCE_TEMPERATURE) ? r_material_properties[REFERENCE_TEMPERATURE] : rValues.GetElementGeometry().GetValue(REFERENCE_TEMPERATURE);
            friction_angle = AdvCLutils::GetPropertyFromTemperatureTable(FRICTION_ANGLE, rValues, ref_temperature);
            cohesion = AdvCLutils::GetPropertyFromTemperatureTable(COHESION, rValues, ref_temperature);
        }
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

        const auto &r_geom = rValues.GetElementGeometry();
        const auto &r_N = rValues.GetShapeFunctionsValues();
        const auto &r_process_info = rValues.GetProcessInfo();

        const double fracture_energy = r_material_properties.GetValue(FRACTURE_ENERGY, r_geom, r_N, r_process_info);
        const double young_modulus   = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);

        double equivalent_yield;
        GetInitialUniaxialThreshold(rValues, equivalent_yield);
        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(equivalent_yield, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture Energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Linear)) { // linear
            rAParameter = -std::pow(equivalent_yield, 2) / (2.0 * young_modulus * fracture_energy / CharacteristicLength);
        } else {
            rAParameter = 0.0;
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
        const array_1d<double, VoigtSize>& rStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rDerivativePlasticPotential,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(rStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
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
        const double friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0;

        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        double c1, c3, c2;
        double checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        if (std::abs(checker) < 29.0) { // If it is not the edge
            c1 = std::sin(friction_angle) / 3.0;
            c3 = (std::sqrt(3.0) * std::sin(lode_angle) + std::sin(friction_angle) * std::cos(lode_angle)) /
                (2.0 * J2 * std::cos(3.0 * lode_angle));
            c2 = 0.5 * std::cos(lode_angle)*(1.0 + std::tan(lode_angle) * std::sin(3.0 * lode_angle) +
                std::sin(friction_angle) * (std::tan(3.0 * lode_angle) - std::tan(lode_angle)) / std::sqrt(3.0));
        } else { // smoothing with drucker-prager
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
        return BaseType::Check(rMaterialProperties);
    }

    /**
     * @brief This method returns true if the yield surfacecompares with the tension tield stress
     */
    static bool IsWorkingWithTensionThreshold()
    {
        return BaseType::IsWorkingWithTensionThreshold();
    }

    /**
     * @brief This method returns the scaling factor of the yield surface surfacecompares with the tension tield stress
     */
    static double GetScaleFactorTension(const Properties& rMaterialProperties)
    {
        return BaseType::GetScaleFactorTension();
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
