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

#if !defined(KRATOS_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED)
#define KRATOS_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// Project includes
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"

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
 * @class MohrCoulombPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a plastic potential following the theory of Mohr-Coulomb
 * @details The Mohr–Coulomb failure surface is a cone with a hexagonal cross section in deviatoric stress space
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo
 */
template <SizeType TVoigtSize = 6>
class MohrCoulombPlasticPotential
{
public:
    ///@name Type Definitions
    ///@{

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// The define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// Counted pointer of MohrCoulombPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombPlasticPotential);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    MohrCoulombPlasticPotential()
    {
    }

    /// Copy constructor
    MohrCoulombPlasticPotential(MohrCoulombPlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    MohrCoulombPlasticPotential &operator=(MohrCoulombPlasticPotential const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~MohrCoulombPlasticPotential(){};

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
		const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double dilatancy = r_material_properties[DILATANCY_ANGLE] * Globals::Pi / 180.0;

        ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        ConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        double c1, c3, c2;
		double checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        if (std::abs(checker) < 29.0) {
            c1 = std::sin(dilatancy);
            c3 = (std::sqrt(3.0) * std::sin(lode_angle) + std::sin(dilatancy) * std::cos(lode_angle)) /
                (2.0 * J2 * std::cos(3.0 * lode_angle));
            c2 = 0.5 * std::cos(lode_angle)*(1.0 + std::tan(lode_angle) * std::sin(3.0 * lode_angle) +
                std::sin(dilatancy) * (std::tan(3.0 * lode_angle) - std::tan(lode_angle)) / std::sqrt(3.0));
        } else { // smoothing with drucker-praguer
            c1 = 3.0 * (2.0 * std::sin(dilatancy) / (std::sqrt(3.0) * (3.0 - std::sin(dilatancy))));
            c2 = 1.0;
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
#endif
