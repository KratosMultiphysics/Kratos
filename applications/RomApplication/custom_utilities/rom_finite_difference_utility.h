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
                                                const ProcessInfo& rCurrentProcessInfo){


        KRATOS_TRY;

        // KRATOS_WATCH("RomFiniteDifferenceUtility::CalculateLeftHandSideDOFDerivative")

        KRATOS_WARNING_IF("RomFiniteDifferenceUtility::CalculateLeftHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
            << "The call of this non omp-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;
        
        #pragma omp critical
        {

            Matrix LHS_p_perturbed;
            Matrix LHS_m_perturbed;

            // Positive perturbation
            // // rElement.InitializeSolutionStep(rCurrentProcessInfo);
            // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            rDof.GetSolutionStepValue() += rPertubationMag;
            // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            // // rElement.FinalizeSolutionStep(rCurrentProcessInfo);
            rElement.CalculateLeftHandSide(LHS_p_perturbed, rCurrentProcessInfo);
            
            // Reset perturbation
            rDof.GetSolutionStepValue() -= rPertubationMag;
            // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);

            // Negative perturbation
            // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            rDof.GetSolutionStepValue() -= rPertubationMag;
            // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            rElement.CalculateLeftHandSide(LHS_m_perturbed, rCurrentProcessInfo);
            
            // Reset perturbation
            rDof.GetSolutionStepValue() += rPertubationMag;
            // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);

            //compute derivative of RHS w.r.t. design variable with finite differences
            rOutput = (LHS_p_perturbed - LHS_m_perturbed) / (2.0*rPertubationMag);
            // KRATOS_WATCH(rOutput)

            // Matrix LHS;
            // Matrix LHS_perturbed;

            // // compute LHS before perturbation
            // rElement.CalculateLeftHandSide(LHS, rCurrentProcessInfo);
            
            // // perturb the design variable and update the element internal values (especially necessary for shell elements with corotational frame)
            // // rElement.InitializeSolutionStep(rCurrentProcessInfo);
            // // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            // rDof.GetSolutionStepValue() += rPertubationMag;
            // // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            // // rElement.FinalizeSolutionStep(rCurrentProcessInfo);

            // // compute LHS after perturbation
            // rElement.CalculateLeftHandSide(LHS_perturbed, rCurrentProcessInfo);

            // //compute derivative of RHS w.r.t. design variable with finite differences
            // rOutput = (LHS_perturbed - LHS) / rPertubationMag;
            // // KRATOS_WATCH(rOutput)
            
            // // unperturb the design variable and update the element internal values (especially necessary for shell elements with corotational frame)
            // // rElement.InitializeSolutionStep(rCurrentProcessInfo);
            // // rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            // rDof.GetSolutionStepValue() -= rPertubationMag;
            // // rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            // // rElement.FinalizeSolutionStep(rCurrentProcessInfo);
            
        }
        KRATOS_CATCH("");
    }

}; // class RomFiniteDifferenceUtility.

}  // namespace Kratos.

#endif // ROM_FINITE_DIFFERENCE_UTILITY_H_INCLUDED  defined


