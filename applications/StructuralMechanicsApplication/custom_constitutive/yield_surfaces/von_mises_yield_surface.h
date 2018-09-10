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

#if !defined(KRATOS_VON_MISES_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_VON_MISES_YIELD_SURFACE_H_INCLUDED

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
 * @class VonMisesYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VonMisesYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of VonMisesYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(VonMisesYieldSurface);

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
        const Vector& rStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double I1, J2;
        Vector Deviator = ZeroVector(6);

        ConstitutiveLawUtilities::CalculateI1Invariant(rStressVector, I1);
        ConstitutiveLawUtilities::CalculateJ2Invariant(rStressVector, I1, Deviator, J2);

        rEquivalentStress = std::sqrt(3.0 * J2);
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
        rThreshold = std::abs(sigma_t);
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param AParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameter(
        ConstitutiveLaw::Parameters& rValues,
        double& AParameter,
        const double CharacteristicLength
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const double Gf = r_material_properties[FRACTURE_ENERGY];
        const double E = r_material_properties[YOUNG_MODULUS];
        const double sigma_c = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];

        if (r_material_properties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            AParameter = 1.00 / (Gf * E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
            KRATOS_ERROR_IF(AParameter < 0.0) << "Fracture enerDerivativePlasticPotentialy is too low, increase FRACTURE_ENERGY..." << std::endl;
        } else { // linear
            AParameter = -std::pow(sigma_c, 2) / (2.0 * E * Gf / CharacteristicLength);
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

        const double c1 = 0.0;
        const double c2 = std::sqrt(3.0);
        const double c3 = 0.0;

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

}; // Class VonMisesYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
