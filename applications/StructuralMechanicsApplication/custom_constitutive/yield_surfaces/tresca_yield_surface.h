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

#if !defined(KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED

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
 * @class TrescaYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Tresca theory
 * @details The Tresca yield criterion is taken to be the work of Henri Tresca. It is also known as the maximum shear stress theory (MSST) and the Tresca–Guest (TG) criterion. 
 * @see https://en.wikipedia.org/wiki/Yield_surface#Tresca_yield_surface
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TrescaYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of TrescaYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(TrescaYieldSurface);
    
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    TrescaYieldSurface()
    {
    }

    /// Copy constructor
    TrescaYieldSurface(TrescaYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    TrescaYieldSurface &operator=(TrescaYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~TrescaYieldSurface(){};

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
        const Vector& rStressVector,
        const Vector& rStrainVector,
        double& rEqStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double I1, J2, J3, lode_angle;
        Vector deviator = ZeroVector(6);

        ConstitutiveLawUtilities::CalculateI1Invariant(rStressVector, I1);
        ConstitutiveLawUtilities::CalculateJ2Invariant(rStressVector, I1, deviator, J2);
        ConstitutiveLawUtilities::CalculateJ3Invariant(deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, lode_angle);

        rEqStress = 2.0 * std::cos(lode_angle) * std::sqrt(J2);
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

        const double sigma_t = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        rThreshold = std::abs(sigma_t); // TODO Check
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
        const double sigma_c = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double sigma_t = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        const double n = sigma_c / sigma_t;

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            rAParameter = 1.00 / (Gf * n * n * E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
            KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture enerDerivativePlasticPotentialy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else { // linear
            rAParameter = -std::pow(sigma_c, 2) / (2.0 * E * Gf * n * n / CharacteristicLength);
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rStressVector The stress vector
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const Vector& rStressVector,
        const Vector& rDeviator,
        const double J2,
        Vector& rDerivativePlasticPotential,
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
     * @param rStressVector The stress vector
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rFFlux The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateYieldSurfaceDerivative(
        const Vector& rStressVector,
        const Vector& rDeviator,
        const double J2,
        Vector& rFFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        Vector first_vector, second_vector, third_vector;

        ConstitutiveLawUtilities::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities::CalculateSecondVector(rDeviator, J2, second_vector);
        ConstitutiveLawUtilities::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities::CalculateJ3Invariant(rDeviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, lode_angle);

        const double checker = std::abs(lode_angle * 180 / Globals::Pi);

        double c2, c3;
        const double c1 = 0.0;

        if (checker < 29.0) {
            c2 = 2.0 * (std::cos(lode_angle) + std::sin(lode_angle) * std::tan(3.0 * lode_angle));
            c3 = std::sqrt(3.0) * std::sin(lode_angle) / (J2 * std::cos(3.0 * lode_angle));
        } else {
            c2 = std::sqrt(3.0);
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

}; // Class TrescaYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
