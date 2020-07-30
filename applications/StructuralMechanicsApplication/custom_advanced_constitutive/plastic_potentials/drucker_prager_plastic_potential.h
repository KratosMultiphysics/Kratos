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

#if !defined(KRATOS_DRUCKER_PRAGER_PLASTIC_POTENTIAL_H_INCLUDED)
#define KRATOS_DRUCKER_PRAGER_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// Project includes
#include "includes/checks.h"
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
 * @class DruckerPragerPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a plastic potential following the theory of Drucker-Prager
 * @details When the yield and plastic potential surfaces are plotted in principal stress space the resulting surface will be a circular cone for Drucker-Prager. This means that both yield and strength are dependent on intermediate principal stress, sigma_2
 * The plastic potential requires the definition of the following properties:
 * - DILATANCY_ANGLE: The angle of dilation controls an amount of plastic volumetric strain developed during plastic shearing and is assumed constant during plastic yielding. The value of DILATANCY_ANGLE=0 corresponds to the volume preserving deformation while in shear.
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <SizeType TVoigtSize = 6>
class DruckerPragerPlasticPotential
{
  public:
    ///@name Type Definitions
    ///@{

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// The define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// Counted pointer of DruckerPragerPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(DruckerPragerPlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    DruckerPragerPlasticPotential()
    {
    }

    /// Copy constructor
    DruckerPragerPlasticPotential(DruckerPragerPlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    DruckerPragerPlasticPotential &operator=(DruckerPragerPlasticPotential const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~DruckerPragerPlasticPotential(){};

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
     * @param J2 The second invariant of the rDeviator
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
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;

        ConstitutiveLawUtilities<VoigtSize>::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);

        const double dilatancy = r_material_properties[DILATANCY_ANGLE] * Globals::Pi / 180.0;
        const double sin_dil = std::sin(dilatancy);
        const double Root3 = std::sqrt(3.0);

        const double CFL = -Root3 * (3.0 - sin_dil) / (3.0 * sin_dil - 3.0);
        const double c1 = CFL * 2.0 * sin_dil / (Root3 * (3.0 - sin_dil));
        const double c2 = CFL;

        noalias(rGFlux) = c1 * first_vector + c2 * second_vector;
    }

    /**
     * @brief This method defines the check to be performed in the plastic potential
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(DILATANCY_ANGLE);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(DILATANCY_ANGLE)) << "DILATANCY_ANGLE is not a defined value" << std::endl;

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
