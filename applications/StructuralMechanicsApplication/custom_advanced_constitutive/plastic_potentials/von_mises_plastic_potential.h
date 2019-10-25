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

#if !defined(KRATOS_VON_MISES_PLASTIC_POTENTIAL_H_INCLUDED)
#define KRATOS_VON_MISES_PLASTIC_POTENTIAL_H_INCLUDED

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
 * @class VonMisesPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a plastic potential following the theory of Von Mises
 * @details If the plastic potential is of vonMises (cylinder) type, on can see that the plastic strain increment tensor is in principle the scaled deviatoric stress tensor, hence principal directions coincide. When the yield and plastic potential surfaces are plotted in principal stress space the resulting surface will be a circular cylinder for Von-Mises. This means that both yield and strength are dependent on intermediate principal stress, sigma_2
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <SizeType TVoigtSize = 6>
class VonMisesPlasticPotential
{
public:
    ///@name Type Definitions
    ///@{

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// The define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// Counted pointer of VonMisesPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(VonMisesPlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    VonMisesPlasticPotential()
    {
    }

    /// Copy constructor
    VonMisesPlasticPotential(VonMisesPlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    VonMisesPlasticPotential &operator=(VonMisesPlasticPotential const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~VonMisesPlasticPotential(){};

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

        ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        const double c1 = 0.0;
        const double c2 = std::sqrt(3.0);
        const double c3 = 0.0;

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
