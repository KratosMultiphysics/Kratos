// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// Project includes
#include "generic_plastic_potential.h"

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
 * @class RankinePlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a plastic potential following the theory of Rankine
 * @author Alejandro Cornejo
 */
template <SizeType TVoigtSize = 6>
class RankinePlasticPotential
{
  public:
    ///@name Type Definitions
    ///@{

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// The define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// Counted pointer of RankinePlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(RankinePlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    RankinePlasticPotential()
    {
    }

    /// Copy constructor
    RankinePlasticPotential(RankinePlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    RankinePlasticPotential &operator=(RankinePlasticPotential const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~RankinePlasticPotential(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This  script  calculates  the derivatives  of the plastic potential
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
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
        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        double c1, c3, c2;
		double checker = std::abs(lode_angle * 180.0 / Globals::Pi);
        const double sqrt_3 = std::sqrt(3.0);

        if (std::abs(checker) < 29.0) { // If it is not the edge
            const double sqrt_J2 = std::sqrt(J2);
            const double square_sin_3_lode = std::pow(std::sin(3.0 * lode_angle), 2);
            const double angle = lode_angle + Globals::Pi / 6.0;
            const double dLode_dJ2 = (3.0 * sqrt_3 * J3) / (4.0 * J2 * J2 * sqrt_J2 * std::sqrt(1.0 - square_sin_3_lode));
            const double dLode_dJ3 = -sqrt_3 / (2.0 * J2 * sqrt_J2 * std::sqrt(1.0 - square_sin_3_lode));
            c1 = 1.0 / 3.0;
            c2 = 2.0 * sqrt_3 / 3.0 * (std::cos(angle) / (2.0 * sqrt_J2) - 2.0 * sqrt_3 * sqrt_J2 / 3.0 * std::sin(angle) * dLode_dJ2) * 2.0 * sqrt_J2;
            c3 = -2.0 * std::sqrt(3.0 * J2) / 3.0 * std::sin(angle) * dLode_dJ3;
        } else { // smoothing with drucker-praguer
            const double friction_angle = rValues.GetMaterialProperties()[FRICTION_ANGLE] * Globals::Pi / 180.0;
            const double sin_phi = std::sin(friction_angle);
            const double CFL = -sqrt_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0);
            c1 = CFL * 2.0 * sin_phi / (sqrt_3 * (3.0 - sin_phi));
            c2 = CFL;
            c3 = 0.0;
        }
        noalias(rGFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
    }

    /**
     * @brief This method defines the check to be performed in the plastic potential
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        return 0;
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

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
