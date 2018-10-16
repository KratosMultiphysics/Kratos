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

#if !defined(KRATOS_DRUCKER_PRAGER_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_DRUCKER_PRAGER_YIELD_SURFACE_H_INCLUDED

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
 * @class DruckerPragerYieldSurface
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
 * @author Alejandro Cornejo & Lucia Barbu
 */
template<class TPlasticPotentialType>
class DruckerPragerYieldSurface
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

    /// Counted pointer of DruckerPragerYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(DruckerPragerYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    DruckerPragerYieldSurface()
    {
    }

    /// Copy constructor
    DruckerPragerYieldSurface(DruckerPragerYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    DruckerPragerYieldSurface &operator=(DruckerPragerYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~DruckerPragerYieldSurface(){};

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
     * @param rEquivalentStress The effective stress or equivalent uniaxial stress is a scalar. It is an invariant value which measures the “intensity” of a 3D stress state.
     */
    static void CalculateEquivalentStress(
        array_1d<double, VoigtSize>& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
        const double sin_phi = std::sin(friction_angle);
        const double root_3 = std::sqrt(3.0);

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("DruckerPragerYieldSurface") << "Friction Angle not defined, assumed equal to 32 " << std::endl;
        }

        double I1, J2;
        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        array_1d<double, VoigtSize> Deviator(VoigtSize, 0.0);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, Deviator, J2);

        if (I1 == 0.0) {
            rEquivalentStress = 0.0;
        } else {
            const double CFL = -root_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0);
            const double TEN0 = 2.0 * I1 * sin_phi / (root_3 * (3.0 - sin_phi)) + std::sqrt(J2);
            rEquivalentStress = std::abs(CFL * TEN0);
        }
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
        const double friction_angle = r_material_properties[FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
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
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const double Gf = r_material_properties[FRACTURE_ENERGY];
        const double E = r_material_properties[YOUNG_MODULUS];
        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double yield_compression = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        const double n = yield_compression / yield_tension;

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (Gf * n * n * E / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else { // linear
            rAParameter = -std::pow(yield_compression, 2) / (2.0 * E * Gf * n * n / CharacteristicLength);
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
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;
        ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        const double c3 = 0.0;

        const double friction_angle = r_material_properties[FRICTION_ANGLE];
        const double sin_phi = std::sin(friction_angle);
        const double Root3 = std::sqrt(3.0);

        const double CFL = -Root3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0);
        const double c1 = CFL * 2.0 * sin_phi / (Root3 * (3.0 - sin_phi));
        const double c2 = CFL;

        noalias(rFFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
    }

    /**
     * @brief This method defines the check to be performed in the yield surface
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(FRICTION_ANGLE);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_TENSION);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_COMPRESSION);
        KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRICTION_ANGLE)) << "FRICTION_ANGLE is not a defined value" << std::endl;
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
     * @brief This method returns true if the yield
	 * surfacecompares with the tension tield stress
     */
    static bool IsWorkingWithTensionThreshold()
    {
        return true;
    }

	/**
     * @brief This method returns the scaling factor of the
     * yield surface
	 * surfacecompares with the tension tield stress
     */
    static double GetScaleFactorTension(const Properties& rMaterialProperties)
    {
        const double friction_angle = rMaterialProperties[FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
        const double sin_phi = std::sin(friction_angle);
        return 1.0 / std::abs((3.0 + sin_phi) / (3.0 * sin_phi - 3.0));
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

}; // Class DruckerPragerYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
