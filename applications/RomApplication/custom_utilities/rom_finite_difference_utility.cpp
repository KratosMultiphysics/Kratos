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


// System includes

// External includes

// Project includes
#include "rom_finite_difference_utility.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
    void CalculateLeftHandSideDOFDerivative(Element& rElement,
                                            Dof<double>& rDof,
                                            const double& rPertubationSize,
                                            Matrix& rOutput,
                                            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_WARNING_IF("RomFiniteDifferenceUtility::CalculateLeftHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
            << "The call of this non omp-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;
        
        #pragma omp critical
        {
            // define working variables
            Matrix LHS;
            Matrix LHS_perturbed;
            Vector dummy;

            // compute LHS before perturbation
            rElement.CalculateLocalSystem(LHS, dummy, rCurrentProcessInfo);
            
            // perturb the design variable
            rDof.GetSolutionStepValue() += rPertubationSize;
            
            // compute LHS after perturbation
            rElement.CalculateLocalSystem(LHS_perturbed, dummy, rCurrentProcessInfo);
            
            //compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (LHS_perturbed - LHS) / rPertubationSize;
            
            // unperturb the design variable
            rDof.GetSolutionStepValue() -= rPertubationSize;

        }
        KRATOS_CATCH("");

    }

}  // namespace Kratos.

