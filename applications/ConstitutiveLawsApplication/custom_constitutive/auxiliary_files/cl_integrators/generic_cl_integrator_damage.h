// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

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
 * @class GenericConstitutiveLawIntegratorDamage
 * @ingroup StructuralMechanicsApplication
 * @brief: This object integrates the predictive stress using the isotropic damage theory by means of
 * linear/exponential softening.
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * The damage integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The fosftening behaviour considered (linear, exponential,etc...)
 * @tparam TYieldSurfaceType The yield surface considered
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TYieldSurfaceType>
class GenericConstitutiveLawIntegratorDamage
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of yield surface
    typedef TYieldSurfaceType YieldSurfaceType;

    /// The define the working dimension size, already defined in the yield surface
    static constexpr SizeType Dimension = YieldSurfaceType::Dimension;

    /// The define the Voigt size, already defined in the yield surface
    static constexpr SizeType VoigtSize = YieldSurfaceType::VoigtSize;

    /// The type of plastic potential
    typedef typename YieldSurfaceType::PlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegratorDamage
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegratorDamage);

    /// Initialization constructor
    GenericConstitutiveLawIntegratorDamage()
    {
    }

    /// Copy constructor
    GenericConstitutiveLawIntegratorDamage(GenericConstitutiveLawIntegratorDamage const &rOther)
    {
    }

    /// Assignment operator
    GenericConstitutiveLawIntegratorDamage &operator=(GenericConstitutiveLawIntegratorDamage const &rOther)
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

    /**
     * @brief This method integrates the predictive stress vector with the CL using linear or exponential softening
     * @param PredictiveStressVector The predictive stress vector
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Damage The internal variable of the damage model
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void IntegrateStressVector(
        array_1d<double, VoigtSize>& rPredictiveStressVector,
        const double UniaxialStress,
        double& rDamage,
        double& rThreshold,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const int softening_type = r_material_properties[SOFTENING_TYPE];
        double damage_parameter;
        TYieldSurfaceType::CalculateDamageParameter(rValues, damage_parameter, CharacteristicLength);

        switch (softening_type)
        {
        case static_cast<int>(SofteningType::Linear):
            CalculateLinearDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        case static_cast<int>(SofteningType::Exponential):
            CalculateExponentialDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        case static_cast<int>(SofteningType::HardeningDamage):
            CalculateHardeningDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        case static_cast<int>(SofteningType::CurveFittingDamage):
            CalculateCurveFittingDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        default:
            KRATOS_ERROR << "SOFTENING_TYPE not defined or wrong..." << softening_type << std::endl;
            break;
        }
        rDamage = (rDamage > 0.99999) ? 0.99999 : rDamage;
        rDamage = (rDamage < 0.0) ? 0.0 : rDamage;
        rPredictiveStressVector *= (1.0 - rDamage);
    }

    /**
     * @brief This computes the damage variable according to exponential softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rDamage The internal variable of the damage model
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateExponentialDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues,
        double& rDamage
        )
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, initial_threshold);
        rDamage = 1.0 - (initial_threshold / UniaxialStress) * std::exp(DamageParameter *
                (1.0 - UniaxialStress / initial_threshold));
    }

    /**
     * @brief This computes the damage variable according to parabolic hardening and exponential
     * softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rDamage The internal variable of the damage model
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateHardeningDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues,
        double& rDamage
        )
    {
        const auto &r_mat_props = rValues.GetMaterialProperties();
        const double max_stress = r_mat_props[MAXIMUM_STRESS];
        const double Gf = r_mat_props[FRACTURE_ENERGY];
        const double E = r_mat_props[YOUNG_MODULUS];
        const bool has_symmetric_yield_stress = r_mat_props.Has(YIELD_STRESS);
        const double yield_compression = has_symmetric_yield_stress ? r_mat_props[YIELD_STRESS] : r_mat_props[YIELD_STRESS_COMPRESSION];
        const double yield_tension = has_symmetric_yield_stress ? r_mat_props[YIELD_STRESS] : r_mat_props[YIELD_STRESS_TENSION];
        const double n = yield_compression / yield_tension;

        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, initial_threshold);

        const double re = max_stress / initial_threshold;
        const double rp = 1.5 * re;
        const double Ad = (rp - re) / re;
        const double Ad_tilda = Ad * (std::pow(rp, 3) - 3.0 * rp + 2.0 / 3.0) / (6.0 * re * std::pow((rp - 1.0), 2));
        const double Hd = 1.0 / (2.0 * (E * Gf * n * n / max_stress / max_stress / CharacteristicLength - 0.5 * rp / re - Ad_tilda));

        const double r = UniaxialStress / initial_threshold;

        if (r <= rp) {
            rDamage = Ad * re / r * std::pow(((r - 1.0) / (rp - 1.0)), 2);
        } else {
            rDamage = 1.0 - re / r + Hd * (1.0 - rp / r);
        }
    }

    /**
     * @brief This computes the damage variable according to linear softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rDamage The internal variable of the damage model
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateLinearDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues,
        double& rDamage
        )
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, initial_threshold);
        rDamage = (1.0 - initial_threshold / UniaxialStress) / (1.0 + DamageParameter);
    }

    /**
     * @brief This computes the damage variable according to a two region curve:
     *          - integrated_stress - strain curve defined by points followed by
     *          - exponential softening curve.
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rDamage The internal variable of the damage model
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateCurveFittingDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues,
        double& rDamage
        )
    {
        const Properties &r_mat_props = rValues.GetMaterialProperties();
        const double fracture_energy = r_mat_props[FRACTURE_ENERGY];
        const double volumetric_fracture_energy = fracture_energy / CharacteristicLength;
        const double yield_stress = r_mat_props[YIELD_STRESS];
        const double E = r_mat_props[YOUNG_MODULUS];
        const Vector& strain_damage_curve = r_mat_props[STRAIN_DAMAGE_CURVE]; //Strain points of the fitting curve
        const Vector& stress_damage_curve = r_mat_props[STRESS_DAMAGE_CURVE]; //Integrated_stress points of the fitting curve
        const SizeType curve_points = strain_damage_curve.size() - 1;

        //Fracture energy required to cover the first region of the curve defined by points
        double volumentric_fracture_energy_first_region = 0.5 * std::pow(yield_stress, 2.0) / E; //Fracture energy corresponding to the elastic regime
        for (IndexType i = 1; i <= curve_points; ++i) {
            volumentric_fracture_energy_first_region += 0.5 * (stress_damage_curve[i-1] + stress_damage_curve[i])
                * (strain_damage_curve[i] - strain_damage_curve[i-1]);
            const double irreversibility_damage_check = (stress_damage_curve[i] - stress_damage_curve[i-1]) / (strain_damage_curve[i] - strain_damage_curve[i-1]);
            KRATOS_ERROR_IF(irreversibility_damage_check > E)<< "The defined S-E curve induces negative damage at region " << i << std::endl;
        }
        KRATOS_ERROR_IF(volumentric_fracture_energy_first_region > volumetric_fracture_energy) << "The Fracture Energy is too low: " << fracture_energy << std::endl;

        const double predictive_stress_end_first_region = strain_damage_curve[curve_points] * E;
        if (UniaxialStress < predictive_stress_end_first_region){ //First region: point-by-point definition with linear interpolation
            for (IndexType i = 1; i <= curve_points; ++i) {
                if (UniaxialStress < strain_damage_curve[i] * E){
                    const double current_integrated_stress = stress_damage_curve[i-1] + (UniaxialStress / E - strain_damage_curve[i-1])
                        * (stress_damage_curve[i] - stress_damage_curve[i-1]) / (strain_damage_curve[i] - strain_damage_curve[i-1]);
                    rDamage = 1.0 - current_integrated_stress / UniaxialStress;
					break;
                }
            }
        } else { //Second region: exponential definition to consume the remaining fracture energy
            const double volumentric_fracture_energy_second_region = volumetric_fracture_energy - volumentric_fracture_energy_first_region;
            rDamage = 1.0 - stress_damage_curve[curve_points] / UniaxialStress * std::exp(stress_damage_curve[curve_points] * (strain_damage_curve[curve_points] * E - UniaxialStress) / (E * volumentric_fracture_energy_second_region));
		}
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThreshold(
        ConstitutiveLaw::Parameters& rValues,
        double& rInitialThreshold
        )
    {
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, rInitialThreshold);
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SOFTENING_TYPE)) << "SOFTENING_TYPE is not a defined value" << std::endl;
        return TYieldSurfaceType::Check(rMaterialProperties);
    }

    /**
     * @brief This method returns the derivative of the yield surface
     * @param rStressVector The stress vector
     * @param rYieldDerivative The derivative of the yield surface
     */
    static void CalculateYieldSurfaceDerivative(
        const array_1d<double, VoigtSize>& rStressVector,
        array_1d<double, VoigtSize>& rYieldDerivative,
        ConstitutiveLaw::Parameters& rValues
    )
    {
        array_1d<double, VoigtSize> deviator = ZeroVector(6);
        double J2;
        const double I1 = rStressVector[0] + rStressVector[1] + rStressVector[2];
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rStressVector, I1, deviator, J2);
        YieldSurfaceType::CalculateYieldSurfaceDerivative(rStressVector, deviator, J2, rYieldDerivative, rValues);
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

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
