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
#include "includes/checks.h"
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
        ConstitutiveLawUtilities<VoigtSize>::CalculateThirdVector(rDeviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ3Invariant(rDeviator, J3);
        ConstitutiveLawUtilities<VoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        const double checker = std::abs(lode_angle * 180.0 / Globals::Pi);

        const double dilatancy = r_material_properties[DILATANCY_ANGLE] * Globals::Pi / 180.0;
        const double sin_dil = std::sin(dilatancy);
        const double cos_dil = std::cos(dilatancy);
        const double sin_theta = std::sin(lode_angle);
        const double cos_theta = std::cos(lode_angle);
        const double cos_3theta = std::cos(3.0 * lode_angle);
        const double tan_theta = std::tan(lode_angle);
        const double tan_3theta = std::tan(3.0 * lode_angle);
        const double Root3 = std::sqrt(3.0);

        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double compr_yield = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double tensi_yield = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        const double n = compr_yield / tensi_yield;

        const double angle_phi = (Globals::Pi * 0.25) + dilatancy * 0.5;
        const double alpha = n / (std::tan(angle_phi) * std::tan(angle_phi));

        const double CFL = 2.0 * std::tan(angle_phi) / cos_dil;

        const double K1 = 0.5 * (1 + alpha) - 0.5 * (1 - alpha) * sin_dil;
        const double K2 = 0.5 * (1 + alpha) - 0.5 * (1 - alpha) / sin_dil;
        const double K3 = 0.5 * (1 + alpha) * sin_dil - 0.5 * (1 - alpha);

        double c1, c2, c3;
        if (std::abs(sin_dil) > tolerance)
            c1 = CFL * K3 / 3.0;
        else
            c1 = 0.0; // check

        if (checker < 29.0) {
            c2 = cos_theta * CFL * (K1 * (1 + tan_theta * tan_3theta) + K2 * sin_dil * (tan_3theta - tan_theta) / Root3);
            c3 = CFL * (K1 * Root3 * sin_theta + K2 * sin_dil * cos_theta) / (2.0 * J2 * cos_3theta);
        } else {
            c3 = 0.0;
            double Aux = 1.0;
            if (std::abs(lode_angle) > tolerance)
                Aux = -1.0;
            c2 = 0.5 * CFL * (K1 * Root3 + Aux * K2 * sin_dil / Root3);
        }

        noalias(rGFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
    }

    /**
     * @brief This method defines the check to be performed in the plastic potential
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(DILATANCY_ANGLE);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_TENSION);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_COMPRESSION);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(DILATANCY_ANGLE)) << "DILATANCY_ANGLE is not a defined value" << std::endl;
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
