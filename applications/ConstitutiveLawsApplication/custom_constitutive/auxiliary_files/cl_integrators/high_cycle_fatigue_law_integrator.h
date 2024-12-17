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
#include "utilities/math_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// #include <numbers>

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
        if (stress_increment_1 > 1.0e-3 && stress_increment_2 < -1.0e-3) {
            rMaximumStress = stress_1;
            rMaxIndicator = true;
        } else if (stress_increment_1 < -1.0e-3 && stress_increment_2 > 1.0e-3) {
            rMinimumStress = stress_1;
            rMinIndicator = true;
        }
    }

    /**
     * @brief This method checks if the global stress state is tension or compression; -1 for a generalized compression state and 1 for a generalized tensile state.
     * @param rStressVector Current predictive stress vector.
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
            const SizeType num_points = stress_damage_curve.size() - 1;

            ultimate_stress = 0.0;
            for (IndexType i = 1; i <= num_points; ++i) {
                ultimate_stress = std::max(ultimate_stress, stress_damage_curve[i-1]);
            }
           
            rUltimateStress = ultimate_stress;
        }
    }

    /**
     * @brief This method computes the stress concentration factor applied to the maximum stress
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param UltimateStress Material ultimate stress.
     * @param rMaterialParameters Material properties.
     * @param rStressConcentrationFactor Stress concentration factor
     */

    static void CalculateStressConcentrationFactor(const double MaxStress,
                                                    const double UltimateStress,
                                                    const Properties& rMaterialParameters,    
                                                    double& rStressConcentrationFactor)
    {       
        const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];
        const double NominalStress = (r_fatigue_coefficients.size() > 9) ? r_fatigue_coefficients[9] * UltimateStress : MaxStress;
        rStressConcentrationFactor = (MaxStress < NominalStress) ? 1.0 : MaxStress / NominalStress;
    }

    /**
     * @brief This method computes the relaxation factor of the residual stresses
     * @param UniaxialStress Equivalent stress of the applied load.
     * @param UniaxialResidualStress Initial equivalent stress.
     * @param InitialThreshold Initial damage threshold.
     * @param rRelaxationFactor Relaxation factor of the residual stresses
     * @param AdvanceStrategyApplied Bool that indicates if the AITS is active
     */

    static void CalculateRelaxationFactor(const double SuperposedUniaxialStress,
                                            const double InitialThreshold,
                                            double& rRelaxationFactor)
    {       
        rRelaxationFactor = (rRelaxationFactor > 1.0) ? 1.0 : (-0.05276 * (SuperposedUniaxialStress / InitialThreshold) + 0.487);
        // rRelaxationFactor = (InitialThreshold - UniaxialStress) / UniaxialResidualStress;  // Relaxation factor proposed by Gustafsson et al. in High cycle fatigue life estimation of punched and trimmed specimens considering residual stress and surface roughness
    }

    /**
     * @brief This method computes the diretional factor at the critical plane to fracture and uniaxial residual stress at this plane
     * @param StressVector previous predictive stress vector.
     * @param DirectionalUniaxialResidualStress uniaxial resitual stress at the critical plane.
     * @param rDirectionalFactor critical plane directional factor.
     */

    static void CalculateDirectionalFactor(const Vector& StressVector,
                                                    const Vector& ResidualStressVector,
                                                    double& rDirectionalUniaxialResidualStress,
                                                    double& rDirectionalFactor)
    {       
        array_1d<double,3> r_n;
        array_1d<double,3> aux_n;

        Matrix stress_tensor(3, 3);
        noalias(stress_tensor) = MathUtils<double>::StressVectorToTensor(StressVector);

        array_1d<double,3> traction_vector;
        array_1d<double,3> tau_vector;
        double sigma_n;
        double tau;
        double directional_uniaxial_stress;
    
        double max_directional_uniaxial_stress = 0.0;

        constexpr SizeType num_points = 2000;
        double golden_angle = M_PI * (3.0 - std::sqrt(5.0));

        for (IndexType i = 0; i < num_points; ++i) {
            // Generate normalized direction vector using Fibonacci sampling
            double z = 1.0 - 2.0 * i / (num_points - 1);  // Map i to [-1, 1]
            double r = std::sqrt(1.0 - std::pow(z, 2.0));
            double theta = i * golden_angle;
            aux_n[0] = r * std::cos(theta);
            aux_n[1] = r * std::sin(theta);
            aux_n[2] = z;       

            // Compute traction vector
            array_1d<double, 3> traction_vector = prod(stress_tensor, aux_n);
            sigma_n = MathUtils<double>::Dot3(aux_n,traction_vector);

            // Compute shear stress           
            noalias(tau_vector) = traction_vector - sigma_n * aux_n;      
            tau = MathUtils<double>::Norm3(tau_vector);

            // Compute von Mises stress at the plane oriented by aux_n normal vector
            directional_uniaxial_stress = std::sqrt(std::pow(sigma_n, 2.0) + 3.0 * std::pow(tau, 2.0));
            
            if (directional_uniaxial_stress > max_directional_uniaxial_stress) {
                max_directional_uniaxial_stress = directional_uniaxial_stress;
                r_n = aux_n;
            }
        }

        const double uniaxial_stress = std::sqrt(0.5 * (std::pow((StressVector[0] - StressVector[1]), 2) + std::pow((StressVector[1] - StressVector[2]), 2) + std::pow((StressVector[2] - StressVector[0]), 2) +
                                                6.0 * (std::pow(StressVector[3], 2) + std::pow(StressVector[4], 2) + std::pow(StressVector[5], 2))));

        rDirectionalFactor = max_directional_uniaxial_stress / uniaxial_stress;

        array_1d<double,3> residual_stress_traction_vector;
        array_1d<double,3> residual_stress_tau_vector;
        double residual_stress_sigma_n;
        double residual_stress_tau;

        Matrix residual_stress_tensor(3, 3);
        noalias(residual_stress_tensor) = MathUtils<double>::StressVectorToTensor(ResidualStressVector);

        noalias(residual_stress_traction_vector) = prod(residual_stress_tensor,r_n);
        residual_stress_sigma_n = MathUtils<double>::Dot3(r_n,residual_stress_traction_vector);

        noalias(residual_stress_tau_vector) = residual_stress_traction_vector - residual_stress_sigma_n * r_n;      
        residual_stress_tau = MathUtils<double>::Norm3(residual_stress_tau_vector);

        rDirectionalUniaxialResidualStress = (residual_stress_sigma_n / std::abs(residual_stress_sigma_n)) * std::sqrt(std::pow(residual_stress_sigma_n, 2.0) + 3.0 * std::pow(residual_stress_tau, 2.0));
    } 

    /**
     * @brief This method computes internal variables (B0, Sth and ALPHAT) of the CL
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param DirectionalUniaxialResidualStress Initial equivalent stress.
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
    static void CalculateFatigueParameters(double MaxStress,
                                            const double DirectionalUniaxialResidualStress,
                                            const double UltimateStress,
                                            double ReversionFactor,
                                            double Threshold,
                                            unsigned int LocalNumberOfCycles,
                                            const Properties& rMaterialParameters,
                                            const ConstitutiveLaw::GeometryType& rElementGeometry,
                                            const double StressConcentrationFactor,
                                            const double DirectionalFactor,
                                            double& rB0,
                                            double& rSth,
                                            double& rAlphat,
                                            double& rN_f,
                                            bool LinearCycleJumpIndicator)
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
        const double FatigueReductionFactorSmoothness = (r_fatigue_coefficients.size() > 7) ? r_fatigue_coefficients[7] : 1.0;
        const double c_stress_concentration = (r_fatigue_coefficients.size() > 8) ? r_fatigue_coefficients[8] : 6.0;
        // const double NominalStress = (r_fatigue_coefficients.size() > 9) ? r_fatigue_coefficients[9] * UltimateStress : MaxStress;

        // Reduction factors applied to the fatigue limit
        
        double nf_trial = LocalNumberOfCycles;
        double reference_max_stress = MaxStress;
        unsigned int max_num_iteration = 50;

        for  (IndexType i = 1; i <= max_num_iteration; ++i){
            
            double h_stress_concentration = (1.0 / (6.0 - c_stress_concentration)) * std::log10(nf_trial / std::pow(10, c_stress_concentration));
            
            if (nf_trial < std::pow(10, c_stress_concentration)){
                h_stress_concentration = 0.0;
            } else if (nf_trial > 1.0e6){
                h_stress_concentration = 1.0;
            }

            MaxStress = reference_max_stress / std::pow(StressConcentrationFactor, (1 - h_stress_concentration));

            const double k_residual_stress = (!DirectionalUniaxialResidualStress) ? 1.0 : (1 - ((DirectionalUniaxialResidualStress + (0.5 + 0.5 * ReversionFactor) * MaxStress * DirectionalFactor) / UltimateStress)) * std::abs((0.5 - 0.5 * ReversionFactor)); // Goodman mean stress correction
            // double const k_residual_stress = 1 - std::pow(((DirectionalUniaxialResidualStress + (0.5 + 0.5 * ReversionFactor) * MaxStress) / UltimateStress), 2.0); // Gerber mean stress correction
            const double k_roughness = (!rElementGeometry.Has(SURFACE_ROUGHNESS)) ? 1.0 : 1 - rElementGeometry.GetValue(MATERIAL_PARAMETER_C1)
                        * std::log10(rElementGeometry.GetValue(SURFACE_ROUGHNESS)) * std::log10((2 * UltimateStress) / rElementGeometry.GetValue(MATERIAL_PARAMETER_C2));
            // const double StressConcentrationFactor = (MaxStress < NominalStress) ? 1.0 : MaxStress / NominalStress;

            if (ReversionFactor <= -1.0) {
                MaxStress *= std::sqrt((1.0 - ReversionFactor) / 2.0); // SWT mean stress correction
                ReversionFactor = - 0.999;
            }

            if (std::abs(ReversionFactor) < 1.0) {
                rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 * (ReversionFactor)), STHR1);
                rAlphat = ALFAF + (0.5 + 0.5 * (ReversionFactor)) * AUXR1;

            } else {
                rSth = Se + (UltimateStress - Se) * std::pow((0.5 + 0.5 / (ReversionFactor)), STHR2);
                rAlphat = ALFAF - (0.5 + 0.5 / (ReversionFactor)) * AUXR2;
            }

            rSth *= (k_residual_stress * k_roughness);
            // MaxStress /= std::pow(StressConcentrationFactor, (1 - h_stress_concentration));

            const double square_betaf = std::pow(BETAF, 2.0);
            double nf_tolerance = 0.1;
            bool is_converged = false;

            if (MaxStress > rSth) {
                if(std::abs(ReversionFactor) < 1.0){
                    rN_f = std::pow(10.0,std::pow(-std::log((MaxStress - rSth) / (UltimateStress - rSth)) / rAlphat,(1.0 / BETAF)));
                    rB0 = -(std::log(MaxStress / UltimateStress) / std::pow((std::log10(rN_f)), FatigueReductionFactorSmoothness * square_betaf));
                    
                    if (std::abs(rN_f - nf_trial) < nf_tolerance){
                        is_converged = true;
                    }
                    
                    const int softening_type = rMaterialParameters[SOFTENING_TYPE];
                    const int curve_by_points = static_cast<int>(SofteningType::CurveFittingDamage);
                    
                    if (is_converged || LinearCycleJumpIndicator){
                        if (softening_type == curve_by_points) {
                            rN_f = std::pow(rN_f, std::pow(std::log(MaxStress / Threshold) / std::log(MaxStress / UltimateStress), 1.0 / (FatigueReductionFactorSmoothness * square_betaf)));
                            if (MaxStress >= Threshold){
                                rN_f = 1.0;
                            }
                        }                           
                        break;
                    }
                    
                    if (std::isnan(rN_f)) {
                        rN_f = std::numeric_limits<double>::infinity();
                        rB0 = -(std::log(MaxStress / UltimateStress) / std::pow((std::log10(rN_f)), FatigueReductionFactorSmoothness * square_betaf));
                        break;
                    }

                    nf_trial = rN_f;

                } else {
                    break;
                } 
            } else {
                break;
            }
            KRATOS_ERROR_IF(i >= max_num_iteration) << "Maximum number of iterations inside the HCF loop exceeded" << std::endl;
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
