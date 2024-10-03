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
     * @param PreviousStresses Equivalent stresses in the two previous steps.
     * @param rMaximumStress Maximum stress.
     * @param rMinimumStress Minimum stress.
     * @param rFirstMaxIndicator Indicator of the first maximun stress found in current cycle.
     * @param rFirstMinIndicator Indicator of the first minimum stress found in current cycle.
     * @param MaxIndicator Indicator of a maximum in the current cycle.
     * @param MinIndicator Indicator of a minimum in the current cycle.
     */
    static void CalculateMaximumAndMinimumStresses(
        const double CurrentStress,
        const Vector& PreviousStresses,
        double& rMaximumStress,
        double& rMinimumStress,
        bool& rFirstMaxIndicator,
        bool& rFirstMinIndicator,
        bool& MaxIndicator,
        bool& MinIndicator)
    {
        const double stress_1 = PreviousStresses[1];
        const double stress_2 = PreviousStresses[0];
        const double stress_increment_1 = stress_1 - stress_2;
        const double stress_increment_2 = CurrentStress - stress_1;
        
        if (stress_increment_1 > 1.0e-3 && stress_increment_2 < -1.0e-3) {
            if (rFirstMaxIndicator){
                rMaximumStress = stress_1;
                rFirstMaxIndicator = false;
                MaxIndicator = true;
            } else if (stress_1 > rMaximumStress){
                rMaximumStress = stress_1;
                MaxIndicator = true;
            }
        } else if (stress_increment_1 < -1.0e-3 && stress_increment_2 > 1.0e-3) {  
            if (rFirstMinIndicator){
                rMinimumStress = stress_1;
                rFirstMinIndicator = false;
                MinIndicator = true;
            } else if (stress_1 < rMinimumStress){
                rMinimumStress = stress_1;
                MinIndicator = true;
            }
        }
    }

    /**
     * @brief This method checks if the global stress state is tension or compression; -1 for a generalized compression state and 1 for a generalized tensile state.
     * @param rStressVector Current predictive stress tensor.
     */
    static double CalculateTensionCompressionFactor(const Vector& rStressVector)
    {
        array_1d<double,3> principal_stresses;
        AdvancedConstitutiveLawUtilities<6>::CalculatePrincipalStresses(principal_stresses, rStressVector);


        double abs_component = 0.0, average_component = 0.0, sum_abs = 0.0, sum_average = 0.0;
        for (unsigned int i = 0; i < principal_stresses.size(); ++i) {
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
     * @brief This method returns the material ultimate stress.
     * @param rUltimateStress Material ultimate stress.
     * @param MaterialParameters Material properties.
     */

    static void CalculateUltimateStress(double& rUltimateStress,
                                          const Properties& rMaterialParameters)
    {
        double ultimate_stress = rMaterialParameters.Has(YIELD_STRESS) ? rMaterialParameters[YIELD_STRESS] : rMaterialParameters[YIELD_STRESS_TENSION];

        // The calculation is prepared to update the rN_f value when using a softening curve which initiates with hardening.
        // The jump in the advance in time process is done in these cases to the Syield rather to Sult.
        const int softening_type = rMaterialParameters[SOFTENING_TYPE];
        const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);
        if (softening_type == curve_by_points) {
            const Vector& stress_damage_curve = rMaterialParameters[STRESS_DAMAGE_CURVE]; //Integrated_stress points of the fitting curve
            const SizeType curve_points = stress_damage_curve.size() - 1;

            ultimate_stress = 0.0;
            for (IndexType i = 1; i <= curve_points; ++i) {
                ultimate_stress = std::max(ultimate_stress, stress_damage_curve[i-1]);
            }
           
            rUltimateStress = ultimate_stress;
        }
    }

    /**
     * @brief This method computes the relaxation factor of the residual stresses
     * @param NominalUniaxialStress Equivalent stress of the applied load.
     * @param UniaxialResidualStress Initial equivalent stress.
     * @param InitialThreshold Initial damage threshold.
     * @param rRelaxationFactor Relaxation factor of the residual stresses
     * @param AdvanceStrategyApplied Bool that indicates if the AITS is active
     */

    static void CalculateRelaxationFactor(const double NominalUniaxialStress,
                                          const double UniaxialResidualStress,    
                                           const double InitialThreshold,
                                           double& rRelaxationFactor)
    {       
        rRelaxationFactor = -0.05276 * ((NominalUniaxialStress + UniaxialResidualStress) / InitialThreshold) + 0.487;
    }

    /**
     * @brief This method computes internal variables (B0, Sth and ALPHAT) of the CL
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param UniaxialResidualStress Initial equivalent stress.
     * @param UltimateStress Material ultimate stress.
     * @param ReversionFactor Ratio between the minimum and maximum signed equivalent stresses for the current load cycle.
     * @param Threshold Damage threshold
     * @param rMaterialParameters Material properties.
     * @param rElementGeometry Element geometry.
     * @param rB0 Internal variable of the fatigue model.
     * @param rSth Endurance limit of the fatigue model.
     * @param rAlphat Internal variable of the fatigue model.
     * @param rN_f Number of cycles to failure.
     */
    static void CalculateFatigueParameters(const double MaxStress,
                                            const double UniaxialResidualStress,
                                            const double UltimateStress,
                                            const double ReversionFactor,
                                            double Threshold,
                                            unsigned int LocalNumberOfCycles,
                                            const Properties& rMaterialParameters,
                                            const ConstitutiveLaw::GeometryType& rElementGeometry,
                                            double& rB0,
                                            double& rSth,
                                            double& rAlphat,
                                            double& rN_f)
	{

        const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];
        const double c1_stress_concentration = r_fatigue_coefficients[8];
        const double c2_stress_concentration = r_fatigue_coefficients[9];
        const double c1_roughness = r_fatigue_coefficients[10];
        const double c2_roughness = r_fatigue_coefficients[11];

        // Reduction factors applied to the fatigue limit
        double const k_residual_stress = 1 - (UniaxialResidualStress / UltimateStress); // Goodman mean stress correction
        // double const k_residual_stress = 1 - std::pow((UniaxialResidualStress / UltimateStress), 2.0); // Gerber mean stress correction
        // double const k_stress_concentration = (!rElementGeometry.Has(HOLE_DIAMETER)) ? 1.0 : 2 + 0.284 * (1 - rElementGeometry.GetValue(HOLE_DIAMETER) / rElementGeometry.GetValue(SPECIMEN_WIDTH)) - 
        //              0.600 * std::pow(1 - rElementGeometry.GetValue(HOLE_DIAMETER) / rElementGeometry.GetValue(SPECIMEN_WIDTH), 2.0) + 1.32 * std::pow(1 - rElementGeometry.GetValue(HOLE_DIAMETER) / rElementGeometry.GetValue(SPECIMEN_WIDTH), 3.0); // Stress concentration factor 
        double const k_stress_concentration = r_fatigue_coefficients[12]; // Manual stress concentration factor
        double const h_stress_concentration = (LocalNumberOfCycles < std::pow(10, c2_stress_concentration)) ? 0.0 : (1.0 / c1_stress_concentration) * std::log10(LocalNumberOfCycles / std::pow(10, c2_stress_concentration));
        double const k_fatigue = (std::pow((1.0 / k_stress_concentration), (h_stress_concentration - 1.0)) > 1.0) ? 1.0 : std::pow((1.0 / k_stress_concentration), (h_stress_concentration - 1.0));
        double const k_roughness = (!rElementGeometry.Has(SURFACE_ROUGHNESS)) ? 1.0 : 1 - rElementGeometry.GetValue(MATERIAL_PARAMETER_C1)
                     * std::log10(rElementGeometry.GetValue(SURFACE_ROUGHNESS)) * std::log10((2 * UltimateStress) / rElementGeometry.GetValue(MATERIAL_PARAMETER_C2));
        double const h_roughness = (LocalNumberOfCycles < std::pow(10, c2_roughness)) ? 0.0 : (1.0 / c1_roughness) * std::log10(LocalNumberOfCycles / std::pow(10, c2_roughness));
                
        //These variables have been defined following the model described by S. Oller et al. in A continuum mechanics model for mechanical fatigue analysis (2005), equation 13 on page 184.
        const double Se = k_residual_stress * k_fatigue * std::pow(k_roughness, h_roughness) * (r_fatigue_coefficients[0] * UltimateStress);
        // const double Se = r_fatigue_coefficients[0] * UltimateStress;
        const double STHR1 = r_fatigue_coefficients[1];
        const double STHR2 = r_fatigue_coefficients[2];
        const double ALFAF = r_fatigue_coefficients[3];
        const double BETAF = r_fatigue_coefficients[4];
        const double AUXR1 = r_fatigue_coefficients[5];
        const double AUXR2 = r_fatigue_coefficients[6];
        const double FatigueReductionFactorSmoothness = r_fatigue_coefficients[7];

        if (std::abs(ReversionFactor) < 1.0) {
            rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 * (ReversionFactor)), STHR1);
			rAlphat = ALFAF + (0.5 + 0.5 * (ReversionFactor)) * AUXR1;
        } else {
            rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 / (ReversionFactor)), STHR2);
			rAlphat = ALFAF - (0.5 + 0.5 / (ReversionFactor)) * AUXR2;
        }

        const double square_betaf = std::pow(BETAF, 2.0);

        if (MaxStress > rSth) {
          if(std::abs(ReversionFactor) < 1.0){
                rN_f = std::pow(10.0,std::pow(-std::log((MaxStress - rSth) / (UltimateStress - rSth)) / rAlphat,(1.0 / BETAF)));
                rB0 = -(std::log(MaxStress / UltimateStress) / std::pow((std::log10(rN_f)), FatigueReductionFactorSmoothness * square_betaf));

                const int softening_type = rMaterialParameters[SOFTENING_TYPE];
                const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);
                
                if (softening_type == curve_by_points) {
                    rN_f = std::pow(rN_f, std::pow(std::log(MaxStress / Threshold) / std::log(MaxStress / UltimateStress), 1.0 / (FatigueReductionFactorSmoothness * square_betaf)));
                    if (MaxStress >= Threshold){
                        rN_f = 1.0;
                    }
                }
                 
                if (std::isnan(rN_f)) {
                    rN_f = std::numeric_limits<double>::infinity();
                    rB0 = -(std::log(MaxStress / UltimateStress) / std::pow((std::log10(rN_f)), FatigueReductionFactorSmoothness * square_betaf));

                    // const double equivalent_max_stress = (0.5 * MaxStress * (1.0 - ReversionFactor)) / (1.0 - (std::pow(((MaxStress * (1.0 + ReversionFactor)) / (2.0 * UltimateStress)), 2.0))); // Gerber mean stress correction
                    // const double reference_reversion_factor = - 0.999;
        
                    // rSth = Se;
                    // rAlphat = ALFAF + (0.5 + 0.5 * (reference_reversion_factor)) * AUXR1;
                    
                    // rN_f = std::pow(10.0,std::pow(-std::log((equivalent_max_stress - rSth) / (UltimateStress - rSth)) / rAlphat,(1.0 / BETAF)));
                    // rB0 = -(std::log(equivalent_max_stress / UltimateStress) / std::pow((std::log10(rN_f)), FatigueReductionFactorSmoothness * square_betaf));

                    // const int softening_type = rMaterialParameters[SOFTENING_TYPE];
                    // const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);

                    // if (softening_type == curve_by_points) {
                    //     rN_f = std::pow(rN_f, std::pow(std::log(equivalent_max_stress / Threshold) / std::log(equivalent_max_stress / UltimateStress), 1.0 / (FatigueReductionFactorSmoothness * square_betaf)));
                    //     if (equivalent_max_stress >= Threshold){
                    //         rN_f = 1.0;
                    //     }
                    // }

                }
        } else if (ReversionFactor <= -1.0) {
                    // const double equivalent_max_stress = (UltimateStress * MaxStress * (1.0 - ReversionFactor)) / (2.0 * UltimateStress - MaxStress * (1.0 + ReversionFactor)); // Goodman mean stress correction
                    const double equivalent_max_stress =  MaxStress * std::sqrt((1.0 - ReversionFactor) / 2.0); // SWT mean stress correction
                    if (equivalent_max_stress > Se) {
                        const double reference_reversion_factor = - 0.999;
           
                        rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 * (ReversionFactor)), STHR1);
                        rAlphat = ALFAF + (0.5 + 0.5 * (reference_reversion_factor)) * AUXR1;
                        
                        rN_f = std::pow(10.0,std::pow(-std::log((equivalent_max_stress - rSth) / (UltimateStress - rSth)) / rAlphat,(1.0 / BETAF)));
                        rB0 = -(std::log(equivalent_max_stress / UltimateStress) / std::pow((std::log10(rN_f)), FatigueReductionFactorSmoothness * square_betaf));

                        const int softening_type = rMaterialParameters[SOFTENING_TYPE];
                        const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);

                        if (softening_type == curve_by_points) {
                            rN_f = std::pow(rN_f, std::pow(std::log(equivalent_max_stress / Threshold) / std::log(equivalent_max_stress / UltimateStress), 1.0 / (FatigueReductionFactorSmoothness * square_betaf)));
                            if (equivalent_max_stress >= Threshold){
                                rN_f = 1.0;
                            }
                        }
                    }
                }
        }
    }

    /**
     * @brief This method computes the reduction factor and the wohler stress (SN curve)
     * @param MaterialParameters Material properties.
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param LocalNumberOfCycles Number of cycles in the current load.
     * @param B0 Internal variable of the fatigue model.
     * @param Sth Endurance limit of the fatigue model.
     * @param rFatigueReductionFactor Reduction factor from the previous step to be reevaluated.
     */
    static void CalculateFatigueReductionFactor(const Properties& rMaterialParameters,
                                                                const double MaxStress,
                                                                unsigned int LocalNumberOfCycles,
                                                                const double B0,
                                                                const double Sth,
                                                                double& rFatigueReductionFactor)
	{
        
        const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];
        const double BETAF = r_fatigue_coefficients[4];
        const double FatigueReductionFactorSmoothness = r_fatigue_coefficients[7];

        if (MaxStress > Sth) {
            rFatigueReductionFactor = std::min(rFatigueReductionFactor, std::exp(-B0 * std::pow(std::log10(static_cast<double>(LocalNumberOfCycles)), FatigueReductionFactorSmoothness * (BETAF * BETAF))));
            // }
            rFatigueReductionFactor = (rFatigueReductionFactor < 0.01) ? 0.01 : rFatigueReductionFactor;
        } 
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
