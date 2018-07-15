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

#if !defined(KRATOS_SIMO_JU_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_SIMO_JU_YIELD_SURFACE_H_INCLUDED

// System includes

// Project includes
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
 * @class SimoJuYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SimoJuYieldSurface
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of SimoJuYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(SimoJuYieldSurface);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    SimoJuYieldSurface()
    {
    }

    /// Copy constructor
    SimoJuYieldSurface(SimoJuYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    SimoJuYieldSurface &operator=(SimoJuYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~SimoJuYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param StressVector The stress vector 
     * @param StrainVector The StrainVector vector
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStress(
        const Vector &StressVector,
        const Vector &StrainVector,
        double &rEqStress,
        const Properties &rMaterialProperties)
    { // It compares with fc / sqrt(E)
        Vector PrincipalStressVector;
        ConstitutiveLawUtilities::CalculatePrincipalStresses(PrincipalStressVector, StressVector);

        double sigma_tension, sigma_compression, n;
        sigma_tension = rMaterialProperties[YIELD_STRESS_TENSION];
        sigma_compression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        n = std::abs(sigma_compression / sigma_tension);

        double SumA = 0.0, SumB = 0.0, SumC = 0.0, ere0, ere1;
        for (std::size_t cont = 0; cont < 2; cont++)
        {
            SumA += std::abs(PrincipalStressVector[cont]);
            SumB += 0.5 * (PrincipalStressVector[cont] + std::abs(PrincipalStressVector[cont]));
            SumC += 0.5 * (-PrincipalStressVector[cont] + std::abs(PrincipalStressVector[cont]));
        }
        ere0 = SumB / SumA;
        ere1 = SumC / SumA;

        double auxf = 0.0;
        for (std::size_t cont = 0; cont < 6; cont++)
        {
            auxf += StrainVector[cont] * StressVector[cont]; // E:S
        }
        rEqStress = std::sqrt(auxf);
        rEqStress *= (ere0 * n + ere1);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rMaterialProperties The material properties
     */
    static void GetInitialUniaxialThreshold(const Properties &rMaterialProperties, double &rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_COMPRESSION] / std::sqrt(rMaterialProperties[YOUNG_MODULUS]));
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param AParameter The damage parameter
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameter(
        const Properties &rMaterialProperties,
        double &AParameter,
        const double CharacteristicLength)
    {
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];
        const double sigma_compression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double sigma_tension = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = sigma_compression / sigma_tension;

        if (rMaterialProperties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential))
        {
            AParameter = 1.0 / (Gf * n * n / (CharacteristicLength * std::pow(sigma_compression, 2)) - 0.5);
            KRATOS_ERROR_IF(AParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        }
        else
        { // linear
            AParameter = -std::pow(sigma_compression, 2) / (2.0 * Gf * n * n / CharacteristicLength);
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param StressVector The stress vector 
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator 
     * @param rg The derivative of the plastic potential
     * @param rMaterialProperties The material properties
     */
    static void CalculatePlasticPotentialDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &rg,
        const Properties &rMaterialProperties)
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, rg, rMaterialProperties);
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param StressVector The stress vector 
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator 
     * @param rFFlux The derivative of the yield surface
     * @param rMaterialProperties The material properties
     */
    static void CalculateYieldSurfaceDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &rFFlux,
        const Properties &rMaterialProperties)
    {
        KRATOS_ERROR << "Yield surface derivative not defined for SimoJu..." << std::endl;
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

}; // Class SimoJuYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
