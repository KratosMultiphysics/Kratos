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

#if !defined(KRATOS_GENERIC_COMPRESSION_CONSTITUTIVE_LAW_INTEGRATOR_DAMAGE_H_INCLUDED)
#define KRATOS_GENERIC_COMPRESSION_CONSTITUTIVE_LAW_INTEGRATOR_DAMAGE_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"

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
 * @class GenericCompressionConstitutiveLawIntegratorDplusDminusDamage
 * @ingroup StructuralMechanicsApplication
 * @brief: This object integrates the predictive stress using the isotropic the d+d- damage theory 
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * The damage integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The fosftening behaviour considered (linear, exponential,etc...)
 * @tparam TYieldSurfaceType The yield surface considered
 * @author Alejandro Cornejo 
 */
template <class TYieldSurfaceType>
class GenericCompressionConstitutiveLawIntegratorDplusDminusDamage
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

    /// Counted pointer of GenericCompressionConstitutiveLawIntegratorDplusDminusDamage
    KRATOS_CLASS_POINTER_DEFINITION(GenericCompressionConstitutiveLawIntegratorDplusDminusDamage);

    /// Initialization constructor
    GenericCompressionConstitutiveLawIntegratorDplusDminusDamage()
    {
    }

    /// Copy constructor
    GenericCompressionConstitutiveLawIntegratorDplusDminusDamage(GenericCompressionConstitutiveLawIntegratorDplusDminusDamage const &rOther)
    {
    }

    /// Assignment operator
    GenericCompressionConstitutiveLawIntegratorDplusDminusDamage &operator=(GenericCompressionConstitutiveLawIntegratorDplusDminusDamage const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericCompressionConstitutiveLawIntegratorDplusDminusDamage()
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
		const int softening_type = r_material_properties.Has(SOFTENING_TYPE_COMPRESSION) ? r_material_properties[SOFTENING_TYPE_COMPRESSION] : r_material_properties[SOFTENING_TYPE];
        double damage_parameter;
        CalculateDamageParameterCompression(rValues, damage_parameter, CharacteristicLength);

        switch (softening_type)
        {
        case static_cast<int>(SofteningType::Linear):
            CalculateLinearDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        case static_cast<int>(SofteningType::Exponential):
            CalculateExponentialDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        default:
            KRATOS_ERROR << "SOFTENING_TYPE not defined or wrong..." << softening_type << std::endl;
            break;
        }
        rPredictiveStressVector *= (1.0 - rDamage);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateDamageParameterCompression(
        ConstitutiveLaw::Parameters& rValues,
        double& rDamageParameter,
        const double CharacteristicLength
        )
    {
        const double fracture_energy_compression = rValues.GetMaterialProperties()[FRACTURE_ENERGY_COMPRESSION];
        ConstitutiveLaw::Parameters modified_values = rValues;
        auto r_properties = modified_values.GetMaterialProperties();
        r_properties.SetValue(FRACTURE_ENERGY, fracture_energy_compression);
        modified_values.SetMaterialProperties(r_properties);
        TYieldSurfaceType::CalculateDamageParameter(modified_values, rDamageParameter, CharacteristicLength);
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
        // This is done to allow tension-driven yields to work as compression yields
        if (YieldSurfaceType::IsWorkingWithTensionThreshold()) {
            ConstitutiveLaw::Parameters modified_ones = rValues;
            const double yield_compression = modified_ones.GetMaterialProperties()[YIELD_STRESS_COMPRESSION];
            Properties material_props = modified_ones.GetMaterialProperties();
            material_props.SetValue(YIELD_STRESS_TENSION, yield_compression);
            modified_ones.SetMaterialProperties(material_props);
            TYieldSurfaceType::GetInitialUniaxialThreshold(modified_ones, rThreshold);  
        } else {
            TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, rThreshold); 
        }   
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
        GetInitialUniaxialThreshold(rValues, initial_threshold);
        rDamage = 1.0 - (initial_threshold / UniaxialStress) * std::exp(DamageParameter *
                 (1.0 - UniaxialStress / initial_threshold));
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
        GetInitialUniaxialThreshold(rValues, initial_threshold);
        rDamage = (1.0 - initial_threshold / UniaxialStress) / (1.0 + DamageParameter);
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(MAXIMUM_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_TENSION);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_COMPRESSION);
        KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
        KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SOFTENING_TYPE)) << "MAXIMUM_STRESS is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS_TENSION)) << "YIELD_STRESS_TENSION is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION)) << "YIELD_STRESS_COMPRESSION is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRACTURE_ENERGY)) << "FRACTURE_ENERGY is not a defined value" << std::endl;

        return TYieldSurfaceType::Check(rMaterialProperties);
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

};
} // namespace Kratos
#endif
