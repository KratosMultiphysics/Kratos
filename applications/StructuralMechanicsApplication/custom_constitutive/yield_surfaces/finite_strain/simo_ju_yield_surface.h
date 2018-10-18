// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//                   Lucia Barbu
//

#if !defined(KRATOS_FINITE_STRAIN_SIMO_JU_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_FINITE_STRAIN_SIMO_JU_YIELD_SURFACE_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h" // TODO: Move to SMALL STRAIN folder
#include "custom_constitutive/yield_surfaces/finite_strain/generic_yield_surface.h"

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
 * @class FiniteStrainSimoJuYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Simo-Ju theory
 * @details The yield surface requires the definition of the following properties:
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Vicente Mataix Ferrandiz
 * @author Alejandro Cornejo
 * @author Lucia Barbu
 */
template <class TPlasticPotentialType>
class FiniteStrainSimoJuYieldSurface
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// The small strain yield surface
    typedef SimoJuYieldSurface<PlasticPotentialType> SmallStrainYieldSurface;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// The definition of the Voigt array type
    typedef array_1d<double, VoigtSize> BoundedArrayType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// Counted pointer of FiniteStrainSimoJuYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(FiniteStrainSimoJuYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    FiniteStrainSimoJuYieldSurface()
    {
    }

    /// Copy constructor
    FiniteStrainSimoJuYieldSurface(FiniteStrainSimoJuYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    FiniteStrainSimoJuYieldSurface &operator=(FiniteStrainSimoJuYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~FiniteStrainSimoJuYieldSurface(){};

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
     * @param rEquivalentStress The effective stress or equivalent uniaxial stress is a scalar. It is an invariant value which measures the “intensity” of a 3D stress state.
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStress(
        const BoundedArrayType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        SmallStrainYieldSurface::CalculateEquivalentStress(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
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
        SmallStrainYieldSurface::GetInitialUniaxialThreshold(rValues, rThreshold);
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
        SmallStrainYieldSurface::CalculateDamageParameter(rValues, rAParameter, CharacteristicLength);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rPredictiveStressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const BoundedArrayType& rPredictiveStressVector,
        const BoundedArrayType& rDeviator,
        const double J2,
        BoundedArrayType& rDerivativePlasticPotential,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rPredictiveStressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rDerivativeYieldSurface The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateYieldSurfaceDerivative(
        const BoundedArrayType& rPredictiveStressVector,
        const BoundedArrayType& rDeviator,
        const double J2,
        BoundedArrayType& rDerivativeYieldSurface,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        SmallStrainYieldSurface::CalculateYieldSurfaceDerivative(rPredictiveStressVector, rDeviator, J2, rDerivativeYieldSurface, rValues);
    }

    /**
     * @brief This method defines the check to be performed in the yield surface
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        return SmallStrainYieldSurface::Check(rMaterialProperties);
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

}; // Class FiniteStrainSimoJuYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
