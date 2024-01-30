// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// Project includes
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"
#include "thermal_von_mises_yield_surface.h"

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
 * @class TrescaYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a yield surface according to Tresca theory with thermal effects
 * @tparam TPlasticPotentialType The plastic potential considered
 * @author Alejandro Cornejo
 */
template <class TPlasticPotentialType>
class ThermalTrescaYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    using PlasticPotentialType = TPlasticPotentialType;

    using BaseType = TrescaYieldSurface<TPlasticPotentialType>;
    using ThermalVonMisesType = ThermalVonMisesYieldSurface<TPlasticPotentialType>;

    /// The Plastic potential already defines the working simension size
    static constexpr SizeType Dimension = PlasticPotentialType::Dimension;

    /// The Plastic potential already defines the Voigt size
    static constexpr SizeType VoigtSize = PlasticPotentialType::VoigtSize;

    /// Counted pointer of TrescaYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ThermalTrescaYieldSurface);

    /// The machine precision zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ThermalTrescaYieldSurface()
    {
    }

    /// Copy constructor
    ThermalTrescaYieldSurface(ThermalTrescaYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    ThermalTrescaYieldSurface &operator=(ThermalTrescaYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ThermalTrescaYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param rStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     * @param rEquivalentStress The effective stress or equivalent uniaxial stress is a scalar. It is an invariant value which measures the “intensity” of a 3D stress state.
     */
    static void CalculateEquivalentStress(
        const array_1d<double, VoigtSize>& rStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BaseType::CalculateEquivalentStress(rStressVector, rStrainVector, rEquivalentStress, rValues);
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
        ThermalVonMisesType::GetInitialUniaxialThreshold(rValues, rThreshold);
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
        ThermalVonMisesType::CalculateDamageParameter(rValues, rAParameter, CharacteristicLength);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param rStressVector The stress vector
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const array_1d<double, VoigtSize>& rStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rDerivativePlasticPotential,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BaseType::CalculatePlasticPotentialDerivative(rStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
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
        const array_1d<double, VoigtSize>& rStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rFFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        BaseType::CalculateYieldSurfaceDerivative(rStressVector, rDeviator, J2, rFFlux, rValues);
    }

    /**
     * @brief This method defines the check to be performed in the yield surface
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        return BaseType::Check(rMaterialProperties);
    }

    /**
     * @brief This method returns true if the yield surfacecompares with the tension tield stress
     */
    static bool IsWorkingWithTensionThreshold()
    {
        return BaseType::IsWorkingWithTensionThreshold();
    }

    /**
     * @brief This method returns the scaling factor of the yield surface surfacecompares with the tension tield stress
     */
    static double GetScaleFactorTension(const Properties& rMaterialProperties)
    {
        return BaseType::GetScaleFactorTension(rMaterialProperties);
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

}; // Class TrescaYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
