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

#if !defined(KRATOS_FINITE_STRAIN_MODIFIED_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED)
#define KRATOS_FINITE_STRAIN_MODIFIED_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h" // TODO: Move to SMALL STRAIN folder
#include "custom_constitutive/plastic_potentials/finite_strain/generic_plastic_potential.h"

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
 * @class FiniteStrainModifiedMohrCoulombPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a plastic potential following the theory of Mohr-Coulomb (modified)
 * @details Working from the conventional assumption that the strength is related to the difference between major and minor principal stresses results in the Tresca model for effective stress. This gives a cone form of the potential in the principal stress space
 * The plastic potential requires the definition of the following properties:
 * - DILATANCY_ANGLE: The angle of dilation controls an amount of plastic volumetric strain developed during plastic shearing and is assumed constant during plastic yielding. The value of DILATANCY_ANGLE=0 corresponds to the volume preserving deformation while in shear.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @author Vicente Mataix Ferrandiz
 * @author Alejandro Cornejo
 * @author Lucia Barbu
 */
template <SizeType TVoigtSize = 6>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiniteStrainModifiedMohrCoulombPlasticPotential
{
  public:
    ///@name Type Definitions
    ///@{

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// The define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// The small strain plastic potential
    typedef ModifiedMohrCoulombPlasticPotential<VoigtSize> SmallStrainPlasticPotential;

    /// The definition of the Voigt array type
    typedef array_1d<double, VoigtSize> BoundedArrayType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// Counted pointer of FiniteStrainModifiedMohrCoulombPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(FiniteStrainModifiedMohrCoulombPlasticPotential);

    /// The zero tolerance definition
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    FiniteStrainModifiedMohrCoulombPlasticPotential()
    {
    }

    /// Copy constructor
    FiniteStrainModifiedMohrCoulombPlasticPotential(FiniteStrainModifiedMohrCoulombPlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    FiniteStrainModifiedMohrCoulombPlasticPotential &operator=(FiniteStrainModifiedMohrCoulombPlasticPotential const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~FiniteStrainModifiedMohrCoulombPlasticPotential(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This  script  calculates  the derivatives  of the plastic potential according   to   NAYAK-ZIENKIEWICZ  paper International journal for numerical methods in engineering vol 113-135 1972.
     * @details As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
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
        SmallStrainPlasticPotential::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
    }

    /**
     * @brief This method defines the check to be performed in the plastic potential
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        return SmallStrainPlasticPotential::Check(rMaterialProperties);
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

}; // Class FiniteStrainModifiedMohrCoulombPlasticPotential

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
