//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, https://github.com/emiroglu
//

#if !defined(ROM_FINITE_DIFFERENCE_UTILITY_H_INCLUDED )
#define  ROM_FINITE_DIFFERENCE_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

/** \brief RomFiniteDifferenceUtility
 *
 * This class calculates the derivatives of different element quantities (e.g. RHS, LHS, mass-matrix, ...)
 * with respect to a design variable (e.g. nodal-coordinate, property).
 */


class KRATOS_API(ROM_APPLICATION) RomFiniteDifferenceUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(RomFiniteDifferenceUtility);

    static void CalculateLeftHandSideDOFDerivative(Element& rElement,
                                                Dof<double>& rDof,
                                                const double& rPertubationMag,
                                                Matrix& rOutput,
                                                bool finite_difference_type,
                                                const ProcessInfo& rCurrentProcessInfo
                                                ){

        KRATOS_TRY;

        KRATOS_WARNING_IF("RomFiniteDifferenceUtility::CalculateLeftHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
            << "The call of this non omp-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;
        
        #pragma omp critical
        {

            // rElement.InitializeSolutionStep(rCurrentProcessInfo);
            // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            // rDof.GetSolutionStepValue() += rPertubationMag;
            // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            // rElement.FinalizeSolutionStep(rCurrentProcessInfo);

            if (finite_difference_type) {

                // Central differencing
                Matrix LHS_p_perturbed;
                Matrix LHS_m_perturbed;
                
                // Positive perturbation
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                rDof.GetSolutionStepValue() += rPertubationMag;
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                rElement.CalculateLeftHandSide(LHS_p_perturbed, rCurrentProcessInfo);
                
                // Negative perturbation
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                rDof.GetSolutionStepValue() -= 2.0*rPertubationMag;
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                rElement.CalculateLeftHandSide(LHS_m_perturbed, rCurrentProcessInfo);

                // Reset perturbation
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                rDof.GetSolutionStepValue() += rPertubationMag;
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                
                // Derivative of LHS w.r.t.
                rOutput = (LHS_p_perturbed - LHS_m_perturbed) / (2.0*rPertubationMag);

            } else {

                // Forward differencing
                Matrix LHS;
                Matrix LHS_p_perturbed;

                // Compute LHS before perturbation
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                rElement.CalculateLeftHandSide(LHS, rCurrentProcessInfo);
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                
                // Positive perturbation
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                rDof.GetSolutionStepValue() += rPertubationMag;
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                rElement.CalculateLeftHandSide(LHS_p_perturbed, rCurrentProcessInfo);
                
                // Reset perturbation
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                rDof.GetSolutionStepValue() -= rPertubationMag;
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                
                // Derivative of LHS w.r.t. DOF
                rOutput = (LHS_p_perturbed - LHS) / rPertubationMag;

            }
                        
        }
        KRATOS_CATCH("");
    }

}; // class RomFiniteDifferenceUtility.

}  // namespace Kratos.

#endif // ROM_FINITE_DIFFERENCE_UTILITY_H_INCLUDED  defined


