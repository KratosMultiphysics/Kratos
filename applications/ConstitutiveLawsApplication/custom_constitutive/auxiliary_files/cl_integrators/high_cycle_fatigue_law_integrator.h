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
//  Main authors:    Sergio Jimenez/Alejandro Cornejo/Lucia Barbu
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    using SizeType = std::size_t;

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
 * @class HighCycleFatigueLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief: This object computes all the required information for the high cycle fatigue constitutive law.
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * @tparam TVoigtSize Strain size
 * @author Sergio Jimenez, Alejandro Cornejo & Lucia Barbu
 */
template <SizeType TVoigtSize = 6>
class HighCycleFatigueLawIntegrator
{
public:
    ///@name Type Definitions
    ///@{

    /// The machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();
    static constexpr double stress_tolerance = 1.0e-3;
    static constexpr SizeType Dimension = (TVoigtSize == 6) ? 3 : 2;

    /// Counted pointer of HighCycleFatigueLawIntegrator
    KRATOS_CLASS_POINTER_DEFINITION(HighCycleFatigueLawIntegrator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor
    HighCycleFatigueLawIntegrator()
    {
    }

    /// Copy constructor
    HighCycleFatigueLawIntegrator(HighCycleFatigueLawIntegrator const &rOther)
    {
    }

    /// Assignment operator
    HighCycleFatigueLawIntegrator &operator=(HighCycleFatigueLawIntegrator const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~HighCycleFatigueLawIntegrator()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks and saves the previous stress state if it was a maximum or a minimum.
     * @param CurrentStress Equivalent stress in the current step.
     * @param rMaximumStress Maximum stress.
     * @param rMinimumStress Minimum stress.
     * @param PreviousStresses Equivalent stresses in the two previous steps.
     * @param rMaxIndicator Indicator of a maximum in the current cycle.
     * @param rMinIndicator Indicator of a minimum in the current cycle.
     */
    static void CalculateMaximumAndMinimumStresses(
        const double CurrentStress,
        double& rMaximumStress,
        double& rMinimumStress,
        const Vector& PreviousStresses,
        bool& rMaxIndicator,
        bool& rMinIndicator)
    {
        const double stress_1 = PreviousStresses[1];
        const double stress_2 = PreviousStresses[0];
        const double stress_increment_1 = stress_1 - stress_2;
        const double stress_increment_2 = CurrentStress - stress_1;
        if (stress_increment_1 > stress_tolerance && stress_increment_2 < -stress_tolerance) {
            rMaximumStress = stress_1;
            rMaxIndicator = true;
        } else if (stress_increment_1 < -stress_tolerance && stress_increment_2 > stress_tolerance) {
            rMinimumStress = stress_1;
            rMinIndicator = true;
        }
    }

    /**
     * @brief This method checks if the global stress state is tension or compression; -1 for a generalized compression state and 1 for a generalized tensile state.
     * @param StressVector Current predictive stress tensor.
     */
    static double CalculateTensionCompressionFactor(const Vector& rStressVector)
    {
        array_1d<double,Dimension> principal_stresses;
        AdvancedConstitutiveLawUtilities<TVoigtSize>::CalculatePrincipalStresses(principal_stresses, rStressVector);


        double abs_component = 0.0, average_component = 0.0, sum_abs = 0.0, sum_average = 0.0;
        for (IndexType i = 0; i < principal_stresses.size(); ++i) {
            abs_component = std::abs(principal_stresses[i]);
            average_component = 0.5 * (principal_stresses[i] + abs_component);
            sum_average += average_component;
            sum_abs += abs_component;
        }
        const double pre_indicator = sum_average / sum_abs;
        if (pre_indicator < 0.5) {
            return -1.0;
        } else {
            return 1.0;
        }
    }

    /**
     * @brief This method returns de reversion factor
     * @param MaxStress Signed maximum equivalent stress in the current cycle.
     * @param MinStress Signed minimum equivalent stress in the current cycle.
     */
    static double CalculateReversionFactor(const double MaxStress, const double MinStress)
    {
        return MinStress / MaxStress;
    }

    /**
     * @brief This method computes internal variables (B0, Sth and ALPHAT) of the CL
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param ReversionFactor Ratio between the minimum and maximum signed equivalent stresses for the current load cycle.
     * @param MaterialParameters Material properties.
     * @param rB0 Internal variable of the fatigue model.
     * @param rSth Endurance limit of the fatigue model.
     * @param rAlphat Internal variable of the fatigue model.
     * @param rNf Number of cycles that satisfy the condition Smax = Su * fred used for the calculation of B0 and, therefore, fred.
     * @param UltimateStress Material ultimate strength.
     * @param FatigueReductionFactorSmoothness Coefficient controlling the smoothness in fred function.
     */
    static void CalculateFatigueParameters(const double MaxStress,
                                            double ReversionFactor,
                                            const Properties& rMaterialParameters,
                                            double& rB0,
                                            double& rSth,
                                            double& rAlphat,
                                            double& rNf,
                                            const double UltimateStress,
                                            const double FatigueReductionFactorSmoothness)
    {
        const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];

        //These variables have been defined following the model described by S. Oller et al. in A continuum mechanics model for mechanical fatigue analysis (2005), equation 13 on page 184.
        const double Se = r_fatigue_coefficients[0] * UltimateStress;
        const double STHR1 = r_fatigue_coefficients[1];
        const double STHR2 = r_fatigue_coefficients[2];
        const double ALFAF = r_fatigue_coefficients[3];
        const double BETAF = r_fatigue_coefficients[4];
        const double AUXR1 = r_fatigue_coefficients[5];
        const double AUXR2 = r_fatigue_coefficients[6];

        if (std::abs(ReversionFactor) < 1.0) {
            rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 * ReversionFactor), STHR1);
			rAlphat = ALFAF + (0.5 + 0.5 * ReversionFactor) * AUXR1;
        } else {
            rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 / ReversionFactor), STHR2);
			rAlphat = ALFAF - (0.5 + 0.5 / ReversionFactor) * AUXR2;
        }

        const double square_betaf = std::pow(BETAF, 2.0);
        if (MaxStress > rSth && MaxStress <= UltimateStress) {
            rNf = std::pow(10.0,std::pow(-std::log((MaxStress - rSth) / (UltimateStress - rSth))/rAlphat,(1.0/BETAF)));
            rB0 = -(std::log(MaxStress / UltimateStress) / std::pow((std::log10(rNf)), FatigueReductionFactorSmoothness * square_betaf));
        } else {
            rNf = std::numeric_limits<double>::infinity(); // No fatigue at this IP, i.e., proposing infinite jump.

        }
    }

    /**
     * @brief This method computes the number of cycles to trigger visible nonlinearities. This will match with the original definition of Nf if
     * the threshold (K) is equal to the ultimate stress (Sult).
     * @param Nf Number of cycles to reactivate the non-linear process.
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param MaterialParameters Material properties.
     * @param Threshold Updated stress threshold at the integration point.
     * @param rSth Endurance limit of the fatigue model.
     * @param UltimateStress Material ultimate strength.
     * @param FatigueReductionFactorSmoothness Coefficient controlling the smoothness in fred function.
     */
    static double NumberOfCyclesToFailure(double Nf,
                                        const double MaxStress,
                                        const Properties& rMaterialParameters,
                                        double Threshold,
                                        double Sth,
                                        const double UltimateStress,
                                        const double FatigueReductionFactorSmoothness)
	{
        const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];
        double number_of_cycles_to_failure = Nf;
        const double BETAF = r_fatigue_coefficients[4];
        if (UltimateStress - MaxStress > tolerance) {
            if (MaxStress > Sth) {
                const double square_betaf = std::pow(BETAF, 2.0);
                number_of_cycles_to_failure = std::pow(Nf, std::pow(std::log(MaxStress / Threshold) / std::log(MaxStress / UltimateStress), 1.0 / (FatigueReductionFactorSmoothness * square_betaf)));
            }
        }
        return number_of_cycles_to_failure;
    }

    /**
     * @brief This method computes the reduction factor and the wohler stress (SN curve)
     * @param MaterialParameters Material properties.
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param LocalNumberOfCycles Number of cycles in the current load.
     * @param GlobalNumberOfCycles Number of cycles in the whole analysis.
     * @param B0 Internal variable of the fatigue model.
     * @param Sth Endurance limit of the fatigue model.
     * @param Alphat Internal variable of the fatigue model.
     * @param rFatigueReductionFactor Reduction factor from the previous step to be reevaluated.
     * @param rWohlerStress Normalized Wohler stress used to build the life prediction curves (SN curve).
     * @param UltimateStress Material ultimate strength.
     * @param FatigueReductionFactorSmoothness Coefficient controlling the smoothness in fred function.
     */
    static void CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters,
                                                                const double MaxStress,
                                                                unsigned int LocalNumberOfCycles,
                                                                unsigned int GlobalNumberOfCycles,
                                                                const double B0,
                                                                const double Sth,
                                                                const double Alphat,
                                                                double& rFatigueReductionFactor,
                                                                double& rWohlerStress,
                                                                const double UltimateStress,
                                                                const double FatigueReductionFactorSmoothness)
	{
        const double BETAF = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS][4];

        if (GlobalNumberOfCycles > 2){
            rWohlerStress = (Sth + (UltimateStress - Sth) * std::exp(-Alphat * (std::pow(std::log10(static_cast<double>(LocalNumberOfCycles)), BETAF)))) / UltimateStress;
        }
        if (MaxStress > Sth) { //In those cases with no fatigue in course (MaxStress < Sth), tbe fatigue reduction factor does not evolve.
            rFatigueReductionFactor = std::min(rFatigueReductionFactor, std::exp(-B0 * std::pow(std::log10(static_cast<double>(LocalNumberOfCycles)), FatigueReductionFactorSmoothness * (BETAF * BETAF))));
            const double min_fatigue_reduction_factor = rMaterialParameters.Has(MINIMUM_FATIGUE_REDUCTION_FACTOR) ? rMaterialParameters[MINIMUM_FATIGUE_REDUCTION_FACTOR] : rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS][0];
            rFatigueReductionFactor = (rFatigueReductionFactor < min_fatigue_reduction_factor) ? min_fatigue_reduction_factor : rFatigueReductionFactor;
        }
    }

    /**
     * @brief This method computes the ultimate stress of the damage model.
     * @param MaterialParameters Material properties.
     */
    static double UltimateStressDamage(const Properties& rMaterialParameters)
	{
        double ultimate_stress = rMaterialParameters.Has(YIELD_STRESS) ? rMaterialParameters[YIELD_STRESS] : rMaterialParameters[YIELD_STRESS_TENSION];

        const int softening_type = rMaterialParameters[SOFTENING_TYPE];
        const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);
        if (softening_type == curve_by_points) {
            const Vector& stress_damage_curve = rMaterialParameters[STRESS_DAMAGE_CURVE]; //Integrated_stress points of the fitting curve
            const SizeType curve_points = stress_damage_curve.size() - 1;

            ultimate_stress = 0.0;
            for (IndexType i = 0; i <= curve_points; ++i) {
                ultimate_stress = std::max(ultimate_stress, stress_damage_curve[i]);
            }
        }
        return ultimate_stress;
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        const double min_fatigue_reduction_factor = rMaterialProperties.Has(MINIMUM_FATIGUE_REDUCTION_FACTOR) ? rMaterialProperties[MINIMUM_FATIGUE_REDUCTION_FACTOR] : rMaterialProperties[HIGH_CYCLE_FATIGUE_COEFFICIENTS][0];
        KRATOS_ERROR_IF(min_fatigue_reduction_factor <= 0.0 || min_fatigue_reduction_factor >= 1.0) << "Minimum fatigue reduction factor must be in (0,1). Provided: " << min_fatigue_reduction_factor << std::endl;

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

}; // Class HighCycleFatigueLawIntegrator

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
