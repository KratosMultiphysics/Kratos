// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// Project includes
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"

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
 * @class ThermalModifiedMohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Modified Mohr-Coulumb theory
 * @details The Mohr–Coulomb yield (failure) criterion is similar to the Tresca criterion, with additional provisions for materials with different tensile and compressive yield strengths. This model is often used to model concrete, soil or granular materials. This is a modified version of the criteria
 * The yield surface requires the definition of the following properties:
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @see https://en.wikipedia.org/wiki/Mohr%E2%80%93Coulomb_theory
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Alejandro Cornejo
 */
template <class TPlasticPotentialType>
class ThermalModifiedMohrCoulombYieldSurface
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    using PlasticPotentialType = TPlasticPotentialType;

    using BaseType = ModifiedMohrCoulombYieldSurface<TPlasticPotentialType>;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// Counted pointer of ThermalModifiedMohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ThermalModifiedMohrCoulombYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Advanced contitutive laws utilities for the corresponding Voigt size
    using AdvCLutils = AdvancedConstitutiveLawUtilities<VoigtSize>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ThermalModifiedMohrCoulombYieldSurface()
    {
    }

    /// Copy constructor
    ThermalModifiedMohrCoulombYieldSurface(ThermalModifiedMohrCoulombYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    ThermalModifiedMohrCoulombYieldSurface &operator=(ThermalModifiedMohrCoulombYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ThermalModifiedMohrCoulombYieldSurface(){};

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
     * @param rEquivalentStress The effective stress or equivalent uniaxial stress is a scalar. It is an invariant value which measures the “intensity” of a 3D stress state.
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStress(
        const array_1d<double, VoigtSize>& rStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        const auto& r_material_properties = rValues.GetMaterialProperties();

        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double yield_compression = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        double friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0; // In radians!

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("ModifiedMohrCoulombYieldSurface") << "Friction Angle not defined, assumed equal to 32 deg " << std::endl;
        }

        double theta;
        const double R = std::abs(yield_compression / yield_tension); // we use the ref ones since the thermal effects acts equally in tension/compression
        const double Rmorh = std::pow(std::tan((Globals::Pi / 4.0) + friction_angle / 2.0), 2);
        const double alpha_r = R / Rmorh;
        const double sin_phi = std::sin(friction_angle);

        double I1, J2, J3;
        array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);
        AdvCLutils::CalculateI1Invariant(rStressVector, I1);
        AdvCLutils::CalculateJ2Invariant(rStressVector, I1, deviator, J2);
        AdvCLutils::CalculateJ3Invariant(deviator, J3);

        const double K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
        const double K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
        const double K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);

        // Check Modified Mohr-Coulomb criterion
        if (std::abs(I1) < tolerance) {
            rEquivalentStress = 0.0;
        } else {
            AdvCLutils::CalculateLodeAngle(J2, J3, theta);
            rEquivalentStress = (2.0 * std::tan(Globals::Pi * 0.25 + friction_angle * 0.5) / std::cos(friction_angle)) * ((I1 * K3 / 3.0) + std::sqrt(J2) * (K1 * std::cos(theta) - K2 * std::sin(theta) * sin_phi / std::sqrt(3.0)));
        }
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThreshold(ConstitutiveLaw::Parameters& rValues, double& rThreshold)
    {
        const auto& r_material_properties = rValues.GetMaterialProperties();

        double yield_compression;
        if (rValues.IsSetShapeFunctionsValues()) { // This is needed since at Initialize level the N are not set yet...
            yield_compression = r_material_properties.Has(YIELD_STRESS) ? AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues) : AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS_COMPRESSION, rValues);
        } else {
            const double ref_temperature = r_material_properties.Has(REFERENCE_TEMPERATURE) ? r_material_properties[REFERENCE_TEMPERATURE] : rValues.GetElementGeometry().GetValue(REFERENCE_TEMPERATURE);
            yield_compression = r_material_properties.Has(YIELD_STRESS) ? AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, ref_temperature) : AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS_COMPRESSION, rValues, ref_temperature);
        }
        rThreshold = std::abs(yield_compression);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviator
     * @param rg The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const array_1d<double, VoigtSize>& rStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& GFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(rStressVector, rDeviator, J2, GFlux, rValues);
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
        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);

        const double yield_compression = has_symmetric_yield_stress ? r_material_properties.GetValue(YIELD_STRESS, r_geom, r_N, r_process_info) : r_material_properties.GetValue(YIELD_STRESS_COMPRESSION, r_geom, r_N, r_process_info);
        const double yield_tension     = has_symmetric_yield_stress ? r_material_properties.GetValue(YIELD_STRESS, r_geom, r_N, r_process_info) : r_material_properties.GetValue(YIELD_STRESS_TENSION, r_geom, r_N, r_process_info);
        const double n = yield_compression / yield_tension;

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (fracture_energy * n * n * young_modulus / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Linear)) { // linear
            rAParameter = -std::pow(yield_compression, 2) / (2.0 * young_modulus * fracture_energy * n * n / CharacteristicLength);
        } else {
            rAParameter = 0.0;
        }
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rStressVector The stress vector
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviator
     * @param rFFlux The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateYieldSurfaceDerivative(
        const array_1d<double, VoigtSize>& rStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rFFlux,
        ConstitutiveLaw::Parameters& rValues)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;

        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        const double checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        double c1, c2, c3;
        double friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0;

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("ModifiedMohrCoulombYieldSurface") << "Friction Angle not defined, assumed equal to 32 deg " << std::endl;
        }

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
        // We use the ref ones since thermal effects are the same in tension and in compression
        const double n = compr_yield / tens_yield;

        const double angle_phi = (Globals::Pi * 0.25) + friction_angle * 0.5;
        const double alpha = n / (std::tan(angle_phi) * std::tan(angle_phi));

        const double CFL = 2.0 * std::tan(angle_phi) / cons_phi;

        const double K1 = 0.5 * (1.0 + alpha) - 0.5 * (1.0 - alpha) * sin_phi;
        const double K2 = 0.5 * (1.0 + alpha) - 0.5 * (1.0 - alpha) / sin_phi;
        const double K3 = 0.5 * (1.0 + alpha) * sin_phi - 0.5 * (1.0 - alpha);

        if (std::abs(sin_phi) > tolerance)
            c1 = CFL * K3 / 3.0;
        else
            c1 = 0.0; // check

        if (std::abs(checker) < 29.0) { // Lode angle not greater than pi/6
            c2 = cos_theta * CFL * (K1 * (1.0 + tan_theta * tan_3theta) + K2 * sin_phi * (tan_3theta - tan_theta) / Root3);
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
        return BaseType::Check(rMaterialProperties);
    }

	/**
     * @brief This method returns true if the yield
	 * surfacecompares with the tension tield stress
     */
    static bool IsWorkingWithTensionThreshold()
    {
        return BaseType::IsWorkingWithTensionThreshold();
    }

	/**
     * @brief This method returns the scaling factor of the
     * yield surface
	 * surfacecompares with the tension tield stress
     */
    static double GetScaleFactorTension(const Properties& rMaterialProperties)
    {
        return BaseType::GetScaleFactorTension(rMaterialProperties);
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

}; // Class ModifiedMohrCoulombYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
