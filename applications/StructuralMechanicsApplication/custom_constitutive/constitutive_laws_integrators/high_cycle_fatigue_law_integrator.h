// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jiménez/Alejandro Cornejo/Lucia Barbu
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
 * @brief: This object integrates the predictive stress using the isotropic damage theory by means of
 * linear/exponential softening.
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * The damage integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The fosftening behaviour considered (linear, exponential,etc...)
 * @tparam TYieldSurfaceType The yield surface considered
 * @author Alejandro Cornejo & Lucia Barbu
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

    static void CalculateMaximumAndMinimumStresses(
        const double CurrentStress,
        double &MaximumStress, 
        double &MinimumStress, 
        const Vector &rPreviousStresses,
        unsigned int& rNumberOfCycles,
        bool& rCycleCounter
        )
    {
        const double Stress1 = rPreviousStresses[1];
        const double Stress2 = rPreviousStresses[0];
		
        const double StressIncrement1 = Stress1 - Stress2;
        const double StressIncrement2 = CurrentStress - Stress1;
        
        if (StressIncrement1 >= 0.001 && StressIncrement2 <= 0.0 && rCycleCounter == false) {
            MaximumStress = Stress1;
            rNumberOfCycles++;
            rCycleCounter = true;
        } else if (StressIncrement1 <= 0.0 && StressIncrement2 >= 0.001 && rCycleCounter == false) {
            MinimumStress = Stress1;
            rCycleCounter = true;
        }
    }

    static void CalculateTensionCompressionFactor(const Vector& StressVector, double& rFactor) 
    {
        array_1d<double,3> principal_stresses;
        ConstitutiveLawUtilities<6>::CalculatePrincipalStresses(principal_stresses, StressVector);


        double abs_component = 0.0, average_component = 0.0, sum_abs = 0.0, sum_average = 0.0;
        for (unsigned int i = 0; i < principal_stresses.size(); ++i) {
            abs_component = std::abs(principal_stresses[i]);
            average_component = 0.5 * (principal_stresses[i] + abs_component);
            sum_average += average_component;
            sum_abs += abs_component;
        }
        const double pre_indicator = sum_average / sum_abs;
        if (pre_indicator < 0.5) {
            rFactor = -1.0;
        } else {
            rFactor = 1.0;
        }
    }

    static void CalculateReversionFactor(const double MaxStress, const double MinStress, double& rReversionFactor)
    {
        rReversionFactor = MinStress / MaxStress;
    }

	static void CalculateFatigueReductionFactor(const double MaxStress,
                                                const double MinStress,
                                                double& rReversionFactor,
                                                const Properties& rMaterialParameters,
                                                const double NumbreOfCycles,
                                                double& rFatigueReductionFactor,
                                                double& rB0
                                                )
	{
        CalculateReversionFactor(MaxStress, MinStress, rReversionFactor);
        const Vector& r_fatigue_parameters = rMaterialParameters[HIGH_CYCLE_FATIGUE_PARAMETERS];
        const double yield_stress = rMaterialParameters.Has(YIELD_STRESS) ? rMaterialParameters[YIELD_STRESS] : rMaterialParameters[YIELD_STRESS_TENSION];

        const double Se = r_fatigue_parameters[0] * yield_stress;
        const double STHR1 = r_fatigue_parameters[1];
        const double STHR2 = r_fatigue_parameters[2];
        const double ALFAF = r_fatigue_parameters[3];
        const double BETAF = r_fatigue_parameters[4];
        const double AUXR1 = r_fatigue_parameters[5];
        const double AUXR2 = r_fatigue_parameters[6];

        double Sth, alphat;
        if (std::abs(rReversionFactor) < 1.0) {
            Sth = Se + (yield_stress - Se) * std::pow((0.5 + 0.5 * rReversionFactor), STHR1);
			alphat = ALFAF + (0.5 + 0.5 * rReversionFactor) * AUXR1;
        } else {
            Sth = Se + (yield_stress - Se) * std::pow((0.5 + 0.5 / rReversionFactor), STHR2);
			alphat = ALFAF - (0.5 + 0.5 / rReversionFactor) * AUXR2;
        }

		KRATOS_WATCH(Sth)
			KRATOS_WATCH(alphat)
			KRATOS_WATCH(rReversionFactor)
        const double square_betaf = std::pow(BETAF, 2.0);
        if (MaxStress > yield_stress) {
            rFatigueReductionFactor = std::exp(-rB0 * std::pow(std::log10(NumbreOfCycles), square_betaf));
        } else if (MaxStress > Sth) {
            const double N_F = std::pow(10,std::pow(-std::log((MaxStress - Sth) / (yield_stress - Sth))/alphat,(1/BETAF)));
            KRATOS_WATCH(std::pow(-std::log((MaxStress - Sth) / (yield_stress - Sth))/alphat,(1/BETAF)))
            //KRATOS_WATCH(std::log((MaxStress - Sth) / (yield_stress - Sth))/alphat)
            //KRATOS_WATCH(1 / BETAF)
			KRATOS_WATCH(N_F)
            rB0 = -(std::log(MaxStress / yield_stress) / std::pow((std::log10(N_F)), square_betaf));
            KRATOS_WATCH(rB0)
            rFatigueReductionFactor = std::exp(-rB0 * std::pow(std::log10(NumbreOfCycles), square_betaf));
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