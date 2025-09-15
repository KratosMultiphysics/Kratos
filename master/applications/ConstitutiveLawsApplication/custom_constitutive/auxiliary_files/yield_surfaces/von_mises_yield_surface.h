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

#pragma once

// System includes

// Project includes
#include "includes/checks.h"
#include "generic_yield_surface.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    using SizeType = std::size_t;

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
 * @class VonMisesYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Von-Mises theory
 * @details The von Mises yield criterion (also known as the maximum distortion energy criterion) suggests that yielding of a ductile material begins when the second deviatoric stress invariant J2 reaches a critical value. It is part of plasticity theory that applies best to ductile materials, such as some metals. Prior to yield, material response can be assumed to be of a nonlinear elastic, viscoelastic, or linear elastic behavior.
 * The yield surface requires the definition of the following properties:
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @see https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class VonMisesYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    using PlasticPotentialType = TPlasticPotentialType;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    using BoundedVector = array_1d<double, VoigtSize>;

    /// Counted pointer of VonMisesYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(VonMisesYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    VonMisesYieldSurface()
    {
    }

    /// Copy constructor
    VonMisesYieldSurface(VonMisesYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    VonMisesYieldSurface &operator=(VonMisesYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~VonMisesYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param rStressVector The stress vector
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
        rEquivalentStress = ConstitutiveLawUtilities<VoigtSize>::CalculateVonMisesEquivalentStress(rStressVector);
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
        const double yield_tension = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
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

        const double fracture_energy = r_material_properties[FRACTURE_ENERGY];
        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        const double yield_compression = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Linear)) { // linear
            rAParameter = -std::pow(yield_compression, 2) / (2.0 * young_modulus * fracture_energy / CharacteristicLength);
        } else {
            rAParameter = 0.0;
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rStressVector The stress vector
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
        const BoundedVector& rPredictiveStressVector,
        const BoundedVector& rDeviator,
        const double J2,
        BoundedVector& rFFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BoundedVector second_vector;
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        const double c2 = std::sqrt(3.0);

        noalias(rFFlux) = c2 * second_vector;
    }

    /**
     * @brief This method defines the check to be performed in the yield surface
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
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

}; // Class VonMisesYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
