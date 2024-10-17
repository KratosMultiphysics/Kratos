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
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"

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
 * @class ThermalDruckerPragerYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Drucker-Prager theory
 * @details The Drucker–Prager yield criterion is similar to the von Mises yield criterion, with provisions for handling materials with differing tensile and compressive yield strengths. This criterion is most often used for concrete where both normal and shear stresses can determine failure.
 * The yield surface requires the definition of the following properties:
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - FRICTION_ANGLE:  Its definition is derived from the Mohr-Coulomb failure criterion and it is used to describe the friction shear resistance of soils together with the normal effective stress.
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @see https://en.wikipedia.org/wiki/Drucker%E2%80%93Prager_yield_criterion
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Alejandro Cornejo
 */
template<class TPlasticPotentialType>
class ThermalDruckerPragerYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    using PlasticPotentialType = TPlasticPotentialType;

    using BaseType = DruckerPragerYieldSurface<TPlasticPotentialType>;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// Counted pointer of ThermalDruckerPragerYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ThermalDruckerPragerYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Advanced contitutive laws utilities for the corresponding Voigt size
    using AdvCLutils = AdvancedConstitutiveLawUtilities<VoigtSize>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ThermalDruckerPragerYieldSurface()
    {
    }

    /// Copy constructor
    ThermalDruckerPragerYieldSurface(ThermalDruckerPragerYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    ThermalDruckerPragerYieldSurface &operator=(ThermalDruckerPragerYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ThermalDruckerPragerYieldSurface(){};

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
     * @param rEquivalentStress The effective stress or equivalent uniaxial stress is a scalar. It is an invariant value which measures the “intensity” of a 3D stress state.
     */
    static void CalculateEquivalentStress(
        array_1d<double, VoigtSize>& rStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {

        double friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0; // In radians!
        const double sin_phi = std::sin(friction_angle);
        const double root_3 = std::sqrt(3.0);

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("DruckerPragerYieldSurface") << "Friction Angle not defined, assumed equal to 32 " << std::endl;
        }

        double I1, J2;
        AdvCLutils::CalculateI1Invariant(rStressVector, I1);
        array_1d<double, VoigtSize> deviator;
        AdvCLutils::CalculateJ2Invariant(rStressVector, I1, deviator, J2);


        const double CFL = -root_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0);
        const double TEN0 = 2.0 * I1 * sin_phi / (root_3 * (3.0 - sin_phi)) + std::sqrt(J2);
        rEquivalentStress = (CFL * TEN0);
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

        double yield_tension, friction_angle;
        if (rValues.IsSetShapeFunctionsValues()) { // This is needed since at Initialize level the N are not set yet...
            yield_tension = r_material_properties.Has(YIELD_STRESS) ? AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues) : AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS_TENSION, rValues);
            friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0;
        } else {
            const double ref_temperature = r_material_properties.Has(REFERENCE_TEMPERATURE) ? r_material_properties[REFERENCE_TEMPERATURE] : rValues.GetElementGeometry().GetValue(REFERENCE_TEMPERATURE);
            yield_tension = r_material_properties.Has(YIELD_STRESS) ? AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, ref_temperature) : AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS_TENSION, rValues, ref_temperature);
            friction_angle = AdvCLutils::GetPropertyFromTemperatureTable(FRICTION_ANGLE, rValues, ref_temperature) * Globals::Pi / 180.0;
        }
        const double sin_phi = std::sin(friction_angle);
        rThreshold = std::abs(yield_tension * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0));
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
        const auto& r_material_properties = rValues.GetMaterialProperties();

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
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rGFlux The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rGFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, rGFlux, rValues);
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
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);

        const double friction_angle = AdvCLutils::GetMaterialPropertyThroughAccessor(FRICTION_ANGLE, rValues) * Globals::Pi / 180.0;;
        const double sin_phi = std::sin(friction_angle);
        const double Root3 = std::sqrt(3.0);

        const double CFL = -Root3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0);
        const double c1 = CFL * 2.0 * sin_phi / (Root3 * (3.0 - sin_phi));
        const double c2 = CFL;

        noalias(rFFlux) = c1 * first_vector + c2 * second_vector;
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
     * yield surface compares with the tension yield stress
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

}; // Class DruckerPragerYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.