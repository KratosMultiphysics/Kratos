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

#include "nlopt_optimizer.h"

namespace Kratos {

NLOptOptimizer::NLOptOptimizer(nlopt_algorithm algorithm, unsigned int num_variables)
    : m_algorithm(algorithm),
      m_num_variables(num_variables),
      m_opt(nullptr),
      m_objective_function(nullptr),
      m_gradient(nullptr),
      m_relative_tolerance(0.0),
      m_absolute_tolerance(0.0),
      m_max_iterations(0)
{
    InitializeOptimizer();
}

NLOptOptimizer::~NLOptOptimizer()
{
    if (m_opt != nullptr) {
        nlopt_destroy(m_opt);
    }
}

void NLOptOptimizer::SetObjectiveFunction(std::function<double(std::vector<double>&)> objective_function)
{
    m_objective_function = objective_function;
}

void NLOptOptimizer::SetGradient(std::function<void(std::vector<double>&, std::vector<double>&)> gradient)
{
    m_gradient = gradient;
}

void NLOptOptimizer::SetLowerBounds(std::vector<double> lower_bounds)
{
    m_lower_bounds = lower_bounds;
}

void NLOptOptimizer::SetUpperBounds(std::vector<double> upper_bounds)
{
    m_upper_bounds = upper_bounds;
}

void NLOptOptimizer::SetInitialGuess(std::vector<double> initial_guess)
{
    m_initial_guess = initial_guess;
}

void NLOptOptimizer::SetRelativeTolerance(double relative_tolerance)
{
    m_relative_tolerance = relative_tolerance;
}

void NLOptOptimizer::SetAbsoluteTolerance(double absolute_tolerance)
{
    m_absolute_tolerance = absolute_tolerance;
}
void NLOptOptimizer::SetMaxIterations(unsigned int max_iterations)
{
    m_max_iterations = max_iterations;
}

std::vector<double> NLOptOptimizer::Optimize()
{
    std::vector<double> result(m_num_variables);

    // Set up the optimization problem
    nlopt_set_min_objective(m_opt, ObjectiveFunctionWrapper, this);
    nlopt_set_lower_bounds(m_opt, m_lower_bounds.data());
    nlopt_set_upper_bounds(m_opt, m_upper_bounds.data());
    nlopt_set_initial_step(m_opt, m_initial_guess.data());
    nlopt_set_xtol_rel(m_opt, m_relative_tolerance);
    nlopt_set_ftol_abs(m_opt, m_absolute_tolerance);
    nlopt_set_maxeval(m_opt, m_max_iterations);

    if (m_gradient != nullptr) {
        //nlopt_set_gradient(m_opt, GradientWrapper, this);
    }

    // Run the optimization algorithm
    double min_value;
    nlopt_result result_code = nlopt_optimize(m_opt, m_initial_guess.data(), &min_value);

    if (result_code < 0) {
        std::string error_message = "NLOpt optimization failed with error code " + std::to_string(result_code);
        throw std::runtime_error(error_message);
    }

    // Copy the optimal values into the result vector
    // const double* optimal_values = nlopt_get_solution(m_opt);
    // std::copy(optimal_values, optimal_values + m_num_variables, result.begin());

    return result;
}

double NLOptOptimizer::ObjectiveFunctionWrapper(unsigned int n, const double* x, double* grad, void* data)
{
    NLOptOptimizer* optimizer = static_cast<NLOptOptimizer*>(data);
    std::vector<double> x_vec(x, x + n);

    double result = optimizer->m_objective_function(x_vec);

    if (grad != nullptr && optimizer->m_gradient != nullptr) {
        optimizer->m_gradient(x_vec, optimizer->m_initial_guess);
        std::copy(optimizer->m_initial_guess.begin(), optimizer->m_initial_guess.end(), grad);
    }

    return result;
}

void NLOptOptimizer::GradientWrapper(unsigned int n, const double* x, double* grad, void* data)
{
    NLOptOptimizer* optimizer = static_cast<NLOptOptimizer*>(data);
    std::vector<double> x_vec(x, x + n);

    optimizer->m_gradient(x_vec, optimizer->m_initial_guess);
    std::copy(optimizer->m_initial_guess.begin(), optimizer->m_initial_guess.end(), grad);
}

void NLOptOptimizer::InitializeOptimizer()
{
    m_opt = nlopt_create(m_algorithm, m_num_variables);

    if (m_opt == nullptr) {
        throw std::runtime_error("Could not create NLOpt optimizer");
    }
}

}