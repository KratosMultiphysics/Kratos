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
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"

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
 * @class ThermalVonMisesYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Von-Mises theory
 * @details The von Mises yield criterion (also known as the maximum distortion energy criterion) suggests that yielding of a ductile material begins when the second deviatoric stress invariant J2 reaches a critical value. It is part of plasticity theory that applies best to ductile materials, such as some metals. Prior to yield, material response can be assumed to be of a nonlinear elastic, viscoelastic, or linear elastic behavior. Includes thermal effects
 * @see https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Alejandro Cornejo
 */
template <class TPlasticPotentialType>
class ThermalVonMisesYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    using PlasticPotentialType = TPlasticPotentialType;

    using BaseType = VonMisesYieldSurface<TPlasticPotentialType>;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// Advanced contitutive laws utilities for the corresponding Voigt size
    using AdvCLutils = AdvancedConstitutiveLawUtilities<VoigtSize>;

    /// Bounded vector for stresses/strains
    using BoundedVector = array_1d<double, VoigtSize>;

    /// Counted pointer of ThermalVonMisesYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ThermalVonMisesYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ThermalVonMisesYieldSurface()
    {
    }

    /// Copy constructor
    ThermalVonMisesYieldSurface(ThermalVonMisesYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    ThermalVonMisesYieldSurface &operator=(ThermalVonMisesYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ThermalVonMisesYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes sqrt(3*J2)
     * @param rStressVector The stress vector
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStress(
        const BoundedVector& rStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BaseType::CalculateEquivalentStress(rStressVector, rStrainVector, rEquivalentStress, rValues);
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
        const auto& r_props = rValues.GetMaterialProperties();
        double yield_tension;
        if (rValues.IsSetShapeFunctionsValues()) { // This is needed since at Initialize level the N are not set yet...
            yield_tension = r_props.Has(YIELD_STRESS) ? AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues) : AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS_TENSION, rValues);
        } else {
            const double ref_temperature = r_props.Has(REFERENCE_TEMPERATURE) ? r_props[REFERENCE_TEMPERATURE] : rValues.GetElementGeometry().GetValue(REFERENCE_TEMPERATURE);
            yield_tension = r_props.Has(YIELD_STRESS) ? AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, ref_temperature) : AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS_TENSION, rValues, ref_temperature);
        }
        rThreshold = std::abs(yield_tension);
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

        // In here we check the accessor
        const double fracture_energy = r_material_properties.GetValue(FRACTURE_ENERGY, r_geom, r_N, r_process_info);
        const double young_modulus   = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);
        const double yield_compression = r_material_properties.Has(YIELD_STRESS) ? r_material_properties.GetValue(YIELD_STRESS, r_geom, r_N, r_process_info) : r_material_properties.GetValue(YIELD_STRESS_COMPRESSION, r_geom, r_N, r_process_info);

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.0 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Linear)) { // linear
            rAParameter = -std::pow(yield_compression, 2) / (2.0 * young_modulus * fracture_energy / CharacteristicLength);
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
        const BoundedVector& rStressVector,
        const BoundedVector& rDeviator,
        const double J2,
        BoundedVector& rDerivativePlasticPotential,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BaseType::CalculatePlasticPotentialDerivative(rStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rStressVector The stress vector
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rFFlux The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateYieldSurfaceDerivative(
        const BoundedVector& rStressVector,
        const BoundedVector& rDeviator,
        const double J2,
        BoundedVector& rFFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BaseType::CalculateYieldSurfaceDerivative(rStressVector, rDeviator, J2, rFFlux, rValues);
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

}; // Class ThermalVonMisesYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
