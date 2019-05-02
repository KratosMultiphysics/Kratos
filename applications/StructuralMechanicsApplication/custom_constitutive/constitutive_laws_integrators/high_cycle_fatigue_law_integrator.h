// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jimenez/Alejandro Cornejo/Lucia Barbu
//

#if !defined(KRATOS_HIGH_CYCLE_FATIGUE_LAW_INTEGRATOR_H_INCLUDED)
#define KRATOS_HIGH_CYCLE_FATIGUE_LAW_INTEGRATOR_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

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
     * @brief This method integrates checks if the uniaxial stress value of the previous step was a maximum or a minimum 
     * @param CurrentStress The signed uniaxial stress in the current step.
     * @param MaximumStress Identifies if the stress in the previous step was a maximum (signed uniaxial stress) in the calculation. Otherwise it is 0.
     * @param MinimumStress Identifies if the stress in the previous step was a minimum (signed uniaxial stress) in the calculation. Otherwise it is 0.
     * @param PreviousMaximumStress Latest identified maximum (signed uniaxial stress) in the calculation.
     * @param PreviousMinimumStress Latest identified minimum (signed uniaxial stress) in the calculation.
     * @param PreviousStresses Vector containing the stress value (signed uniaxial stresses) from the two previous steps.
     * @param NumberOfCycles Number of cycles for the applied load.
     * @param CycleCounter Boolean variable used on the update operation for the SetMaxStress.
     */
    static void CalculateMaximumAndMinimumStresses(
        const double CurrentStress,
        double &rMaximumStress, 
        double &rMinimumStress,  
        double &rPreviousMaximumStress,
        double &rPreviousMinimumStress,
        const Vector &rPreviousStresses,
        unsigned int& rNumberOfCycles,
        bool& rCycleCounter
        )
    {
        const double stress_1 = rPreviousStresses[1];
        const double stress_2 = rPreviousStresses[0];
        const double stress_increment_1 = stress_1 - stress_2;
        const double stress_increment_2 = CurrentStress - stress_1;

        if (stress_increment_1 >= 0.001 && stress_increment_2 <= 0.0 && !rCycleCounter) {
            rMaximumStress = stress_1;
            rPreviousMaximumStress = rMaximumStress;
            rCycleCounter = true;
        } else if (stress_increment_1 <= -0.001 && stress_increment_2 >= 0.0 && !rCycleCounter) {
            rMinimumStress = stress_1;
            rPreviousMinimumStress = rMinimumStress;
            rCycleCounter = true;

        }
        
        if (rPreviousMaximumStress != rMaximumStress && rPreviousMinimumStress != rMinimumStress && !rCycleCounter) {
            rNumberOfCycles++;
            rPreviousMaximumStress = 0.0;
            rPreviousMinimumStress = 0.0;
        }
    }

    /**
     * @brief This method checks if the global stress state is tension or compression; -1 for a generalized compression state and 1 for a generalized tensile state. 
     * @param StressVector Current predictive stress tensor.
     */
    static double CalculateTensionCompressionFactor(const Vector& rStressVector) 
    {
        array_1d<double,3> principal_stresses;
        ConstitutiveLawUtilities<6>::CalculatePrincipalStresses(principal_stresses, rStressVector);


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
     * @brief This method computes the reduction factor used in the high cycle fatigue model 
     * @param MaxStress Signed maximum stress in the current cycle.
     * @param MinStress Signed minimum stress in the current cycle.
     * @param ReversionFactor Ratio between the minimum and maximum signed equivalent stresses for the current load cycle.
     * @param MaterialParameters Material properties.
     * @param NumbreOfCycles Number of cycles at the current step.
     * @param FatigueReductionFactor Reduction factor from the previous step to be reevaluated.
     * @param B0 Internal variable of the fatigue model used for the calculation of the FatigueReductionFactor.
     */
    static void CalculateFatigueReductionFactor(const double MaxStress,
                                                const double MinStress,
                                                double ReversionFactor,
                                                const Properties& rMaterialParameters,
                                                const int NumberOfCycles,
                                                double& rFatigueReductionFactor,
                                                double& rB0
                                                )
	{
        const Vector& r_fatigue_coefficients = rMaterialParameters[HIGH_CYCLE_FATIGUE_COEFFICIENTS];
        const double yield_stress = rMaterialParameters.Has(YIELD_STRESS) ? rMaterialParameters[YIELD_STRESS] : rMaterialParameters[YIELD_STRESS_TENSION];

        //These variables have been defined following the model described by S. Oller et al. in A continuum mechanics model for mechanical fatigue analysis (2005), equation 13 on page 184.
        const double Se = r_fatigue_coefficients[0] * yield_stress;
        const double STHR1 = r_fatigue_coefficients[1];  
        const double STHR2 = r_fatigue_coefficients[2];
        const double ALFAF = r_fatigue_coefficients[3];
        const double BETAF = r_fatigue_coefficients[4];
        const double AUXR1 = r_fatigue_coefficients[5];
        const double AUXR2 = r_fatigue_coefficients[6];

        double Sth, alphat;
        if (std::abs(ReversionFactor) < 1.0) {
            Sth = Se + (yield_stress - Se) * std::pow((0.5 + 0.5 * ReversionFactor), STHR1);
			alphat = ALFAF + (0.5 + 0.5 * ReversionFactor) * AUXR1;
        } else {
            Sth = Se + (yield_stress - Se) * std::pow((0.5 + 0.5 / ReversionFactor), STHR2);
			alphat = ALFAF - (0.5 + 0.5 / ReversionFactor) * AUXR2;
        }

        const double square_betaf = std::pow(BETAF, 2.0);
        if (MaxStress > yield_stress) {
            rFatigueReductionFactor = std::exp(-rB0 * std::pow(std::log10(static_cast<double>(NumberOfCycles)), square_betaf));
        } else if (MaxStress > Sth) {
            const double N_F = std::pow(10.0,std::pow(-std::log((MaxStress - Sth) / (yield_stress - Sth))/alphat,(1.0/BETAF)));
            rB0 = -(std::log(MaxStress / yield_stress) / std::pow((std::log10(N_F)), square_betaf));
            rFatigueReductionFactor = std::exp(-rB0 * std::pow(std::log10(static_cast<double>(NumberOfCycles)), square_betaf));
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
#endif