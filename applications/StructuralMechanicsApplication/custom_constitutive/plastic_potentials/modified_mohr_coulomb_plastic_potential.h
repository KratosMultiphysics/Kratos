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

#if !defined(KRATOS_MODIFIED_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED)
#define KRATOS_MODIFIED_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"

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
 * @class ModifiedMohrCoulombPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ModifiedMohrCoulombPlasticPotential
{
  public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ModifiedMohrCoulombPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(ModifiedMohrCoulombPlasticPotential);

    double tolerance = 100000000;// std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ModifiedMohrCoulombPlasticPotential()
    {
    }

    /// Copy constructor
    ModifiedMohrCoulombPlasticPotential(ModifiedMohrCoulombPlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    ModifiedMohrCoulombPlasticPotential &operator=(ModifiedMohrCoulombPlasticPotential const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ModifiedMohrCoulombPlasticPotential(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This  script  calculates  the derivatives  of the plastic potential according   to   NAYAK-ZIENKIEWICZ  paper International journal for numerical methods in engineering vol 113-135 1972.
     * @details As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param StressVector The stress vector 
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator 
     * @param rGFlux The derivative of the plastic potential
     * @param rMaterialProperties The material properties
     */
    static void CalculatePlasticPotentialDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &rGFlux,
        const Properties &rMaterialProperties)
    {
        Vector first_vector, second_vector, third_vector;

        ConstitutiveLawUtilities::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities::CalculateSecondVector(Deviator, J2, second_vector);
        ConstitutiveLawUtilities::CalculateThirdVector(Deviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, lode_angle);

        const double Checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        const double dilatancy = rMaterialProperties[DILATANCY_ANGLE] * Globals::Pi / 180.0;
        const double sin_dil = std::sin(dilatancy);
        const double cos_dil = std::cos(dilatancy);
        const double sin_theta = std::sin(lode_angle);
        const double cos_theta = std::cos(lode_angle);
        const double cos_3theta = std::cos(3.0 * lode_angle);
        const double tan_theta = std::tan(lode_angle);
        const double tan_3theta = std::tan(3.0 * lode_angle);
        const double Root3 = std::sqrt(3.0);

        const double compr_yield = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double tensi_yield = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = compr_yield / tensi_yield;

        const double angle_phi = (Globals::Pi * 0.25) + dilatancy * 0.5;
        const double alpha = n / (std::tan(angle_phi) * std::tan(angle_phi));

        const double CFL = 2.0 * std::tan(angle_phi) / cos_dil;

        const double K1 = 0.5 * (1 + alpha) - 0.5 * (1 - alpha) * sin_dil;
        const double K2 = 0.5 * (1 + alpha) - 0.5 * (1 - alpha) / sin_dil;
        const double K3 = 0.5 * (1 + alpha) * sin_dil - 0.5 * (1 - alpha);

        double c1, c2, c3;
        if (std::abs(sin_dil) > 0.01)
            c1 = CFL * K3 / 3.0;
        else
            c1 = 0.0; // check

        if (Checker < 29.0)
        {
            c2 = cos_theta * CFL * (K1 * (1 + tan_theta * tan_3theta) + K2 * sin_dil * (tan_3theta - tan_theta) / Root3);
            c3 = CFL * (K1 * Root3 * sin_theta + K2 * sin_dil * cos_theta) / (2.0 * J2 * cos_3theta);
        }
        else
        {
            c3 = 0.0;
            double Aux = 1.0;
            if (std::abs(lode_angle) > 0.01)
                Aux = -1.0;
            c2 = 0.5 * CFL * (K1 * Root3 + Aux * K2 * sin_dil / Root3);
        }

        noalias(rGFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
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

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
