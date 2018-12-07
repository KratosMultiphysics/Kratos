//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/convergence_accelerator.hpp"
#include "custom_utilities/mvqn_convergence_accelerator.hpp"
#include "custom_utilities/mvqn_recursive_convergence_accelerator.hpp"

namespace Kratos {
	namespace Testing {

        /** 
	     * Auxiliar function to set the system to be solved
	     */
        void ComputeResidual(
            const std::size_t TimeValue,
            Vector& rGuess,
            Vector& rRes) {

            // Set A matrix and b (RHS) vector
            Vector b = ZeroVector(5);
            Matrix A = ZeroMatrix(5,5);
            b(0) = std::cos(0.05*TimeValue*Globals::Pi);
            b(1) = 2.0 * std::cos(0.05*TimeValue*Globals::Pi);
            b(2) = 3.0 * std::cos(0.05*TimeValue*Globals::Pi);
            b(3) = 4.0 * std::cos(0.05*TimeValue*Globals::Pi);
            b(4) = 5.0 * std::cos(0.05*TimeValue*Globals::Pi);
            A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
            A(1,1) = 4.0; A(1,2) = 5.0; A(1,3) = 6.0;
            A(2,2) = 7.0; A(2,3) = 8.0; A(2,4) = 9.0;
            A(3,3) = 10.0; A(4,4) = 11.0;
            A(4,0) = 1.0; A(4,1) = 2.0; A(4,2) = 3.0; A(4,3) = 4.0; A(4,4) = 5.0;

            // Compute the residual vector
            Vector aux = ZeroVector(5);
            UblasSpace<double, Matrix, Vector >::Mult(A, rGuess, aux);
            UblasSpace<double, Matrix, Vector >::Assign(rRes, 1.0, b);
            UblasSpace<double, Matrix, Vector >::UnaliasedAdd(rRes, -1.0, aux);
        }

        template<class TSpace>
        bool SolveProblem(
            typename ConvergenceAccelerator<TSpace>::Pointer pConvAccelerator,
            const double Tol = 1e-9,
            const std::size_t MaxIt = 25,
            const std::size_t EndTime = 10) {

            // Initialize the initial guess and residual vectors
            Vector res = ZeroVector(5);
            Vector guess = ZeroVector(5);
            for (std::size_t i = 0; i < 5; ++i){
                res(i) = 1.0;
                guess(i) = 1.0;
            }

            // Initialize the convergence accelerator
            pConvAccelerator->Initialize();

            // Perform the iteration until convergence
            double res_norm = 0.0;
            for (std::size_t t_val = 0; t_val < EndTime; ++t_val){
                std::cout << "\nStep " << t_val + 1 << " resolution starts..." << std::endl;
                pConvAccelerator->InitializeSolutionStep();
                for (std::size_t it = 0; it < MaxIt; ++it){
                    ComputeResidual(t_val, guess, res);
                    res_norm = TSpace::TwoNorm(res);
                    std::cout << "\tIteration: " << it + 1 << " residual: " << res_norm << std::endl; 
                    if (res_norm < Tol){
                        std::cout << "Convergence achieved in " << it + 1 << " iterations." << std::endl;
                        break;
                    } else {
                        pConvAccelerator->InitializeNonLinearIteration();
                        pConvAccelerator->UpdateSolution(res, guess);
                        pConvAccelerator->FinalizeNonLinearIteration();
                    }
                }
                pConvAccelerator->FinalizeSolutionStep();
            }

            // Return true if the problem has converged
            return (res_norm < Tol) ? true : false;
        }

	    /** 
	     * Checks the MVQN convergence accelerator
	     */
	    KRATOS_TEST_CASE_IN_SUITE(MVQNConvergenceAccelerator, FSIApplicationFastSuite)
		{
            // Set the convergence accelerator pointer
            const double w_0 = 0.825;
            MVQNFullJacobianConvergenceAccelerator<UblasSpace<double, Matrix, Vector > >::Pointer pMVQN = 
                Kratos::make_shared<MVQNFullJacobianConvergenceAccelerator<UblasSpace<double, Matrix, Vector > > >(w_0);

            // Solve the Ax = b problem
            const bool is_converged = SolveProblem<UblasSpace<double, Matrix, Vector > >(pMVQN);

            // Check results
            KRATOS_CHECK(is_converged);
	    }

	    /** 
	     * Checks the recursive MVQN convergence accelerator
	     */
	    KRATOS_TEST_CASE_IN_SUITE(RecursiveMVQNConvergenceAccelerator, FSIApplicationFastSuite)
		{
            // Set the convergence accelerator
            const double w_0 = 0.825;
            const std::size_t buffer_size = 10;
            MVQNRecursiveJacobianConvergenceAccelerator<UblasSpace<double, Matrix, Vector > >::Pointer pRecursiveMVQN =
                Kratos::make_shared<MVQNRecursiveJacobianConvergenceAccelerator<UblasSpace<double, Matrix, Vector > > >(w_0, buffer_size);

            // Solve the Ax = b problem
            const bool is_converged = SolveProblem< UblasSpace<double, Matrix, Vector > >(pRecursiveMVQN);

            // Check results
            KRATOS_CHECK(is_converged);
	    }

    } // namespace Testing
}  // namespace Kratos.
