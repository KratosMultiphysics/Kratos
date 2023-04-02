//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#ifndef NLOPT_OPTIMIZER_H
#define NLOPT_OPTIMIZER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "algorithm_base.h"

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include "custom_external_libraries/nlopt/src/api/nlopt.h"

// ==============================================================================

namespace Kratos {

class KRATOS_API(OPTIMIZATION_APPLICATION) NLOptOptimizer 
{
    public:
        KRATOS_CLASS_POINTER_DEFINITION(NLOptOptimizer);

        NLOptOptimizer(nlopt_algorithm algorithm, unsigned int num_variables);
        virtual ~NLOptOptimizer();
        void SetObjectiveFunction(std::function<double(std::vector<double>&)> objective_function);
        void SetGradient(std::function<void(std::vector<double>&, std::vector<double>&)> gradient);
        void SetLowerBounds(std::vector<double> lower_bounds);
        void SetUpperBounds(std::vector<double> upper_bounds);
        void SetInitialGuess(std::vector<double> initial_guess);
        void SetRelativeTolerance(double relative_tolerance);
        void SetAbsoluteTolerance(double absolute_tolerance);
        void SetMaxIterations(unsigned int max_iterations);

        std::vector<double> Optimize();

        private:
            static double ObjectiveFunctionWrapper(unsigned n, const double* x, double* grad, void* data);
            static void GradientWrapper(unsigned n, const double* x, double* grad, void* data);
            void InitializeOptimizer();

            nlopt_algorithm m_algorithm;
            unsigned int m_num_variables;
            nlopt_opt m_opt;
            std::function<double(std::vector<double>&)> m_objective_function;
            std::function<void(std::vector<double>&, std::vector<double>&)> m_gradient;
            std::vector<double> m_lower_bounds;
            std::vector<double> m_upper_bounds;
            std::vector<double> m_initial_guess;
            double m_relative_tolerance;
            double m_absolute_tolerance;
            unsigned int m_max_iterations;

    }; // Class NLOptOptimizer

} // namespace Kratos

#endif // NLOPT_OPTIMIZER_H
