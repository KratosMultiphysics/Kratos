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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericCompressionConstitutiveLawIntegratorDplusDminusDamage
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
        const double peak_stress_compression = r_material_properties[MAXIMUM_STRESS];

        double initial_threshold;
        this->GetInitialUniaxialThreshold(rValues, initial_threshold);

        const double Ad = (peak_stress_compression - initial_threshold) / initial_threshold;

        if (UniaxialStress <= peak_stress_compression) { // Polinomic path
            rDamage = Ad * (initial_threshold / UniaxialStress) * ((UniaxialStress - initial_threshold) / (peak_stress_compression - initial_threshold));
        } else { // Exponential softening
            const double Gf = r_material_properties[FRACTURE_ENERGY];
            const double E = r_material_properties[YOUNG_MODULUS];
            const double Ad_hat = Ad * (std::pow(peak_stress_compression, 3) - 3.0 * peak_stress_compression * std::pow(initial_threshold, 2) + 2.0 * std::pow(initial_threshold, 3)) /
                                  (6.0 * initial_threshold * std::pow((peak_stress_compression - initial_threshold), 2));
            const double Hd = 0.5 / (E * 100 * Gf / initial_threshold / CharacteristicLength - 0.5 * peak_stress_compression / initial_threshold - Ad_hat);
            rDamage = 1.0 - initial_threshold / UniaxialStress * std::exp(2.0 * Hd * (peak_stress_compression - UniaxialStress) / initial_threshold);
        }
        rPredictiveStressVector *= (1.0 - rDamage);
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
        ConstitutiveLaw::Parameters modified_ones = rValues;

        // This is done to allow tension-driven yields to work as compression yields
        if (TConstLawIntegratorCompressionType::YieldSurfaceType::IsWorkingWithTensionThreshold()) {
            const double yield_compression = modified_ones.GetMaterialProperties()[YIELD_STRESS_COMPRESSION];
            Properties material_props = modified_ones.GetMaterialProperties();
            material_props.SetValue(YIELD_STRESS_TENSION, yield_compression);
            modified_ones.SetMaterialProperties(material_props);
        }
        TConstLawIntegratorCompressionType::GetInitialUniaxialThreshold(modified_ones, initial_threshold_compression);        
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(MAXIMUM_STRESS);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SOFTENING_TYPE)) << "MAXIMUM_STRESS is not a defined value" << std::endl;

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

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const
    {
    }

    void load(Serializer &rSerializer)
    {
    }

    ///@}

};
} // namespace Kratos
#endif