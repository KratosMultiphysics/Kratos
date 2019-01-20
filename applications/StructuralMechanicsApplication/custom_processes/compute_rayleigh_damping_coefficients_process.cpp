// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/linear_solver_factory.h"
#include "utilities/math_utils.h"
#include "custom_processes/compute_rayleigh_damping_coefficients_process.h"
#include "structural_mechanics_application_variables.h"
#include "custom_strategies/custom_schemes/eigensolver_dynamic_scheme.hpp"
#include "custom_strategies/custom_strategies/eigensolver_strategy.hpp"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"

namespace Kratos
{

void ComputeRayleighDampingCoefficientsProcess::Execute()
{
    KRATOS_TRY

    // Getting process info
    ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

    // Getting echo level
    const int echo_level = mParameters["echo_level"].GetInt();

    // If not computed eigen values already
    if (!r_current_process_info.Has(EIGENVALUE_VECTOR)) {
        // Getting Feast linear solver
        typedef TUblasSparseSpace<double> SparseSpaceType;
        typedef TUblasDenseSpace<double> LocalSparseSpaceType;
        typedef LinearSolver<SparseSpaceType,LocalSparseSpaceType> LinearSolverType;
        typedef LinearSolverFactory<SparseSpaceType,  LocalSparseSpaceType> LinearSolverFactoryType;

        LinearSolverType::Pointer p_linear_solver = LinearSolverFactoryType().Create(mParameters["eigen_system_settings"]);

        // Getting B&S
        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolverWithConstraints<SparseSpaceType, LocalSparseSpaceType, LinearSolverType>>(p_linear_solver);

        // Getting scheme
        auto p_scheme = Kratos::make_shared<EigensolverDynamicScheme<SparseSpaceType, LocalSparseSpaceType>>();

        // Getting strategy
        auto p_compute_eigen_values_strategy = Kratos::make_shared<EigensolverStrategy<SparseSpaceType, LocalSparseSpaceType, LinearSolverType>>(mrModelPart, p_scheme, p_builder_and_solver);

        // Compute the eigen values
        KRATOS_INFO_IF("ComputeRayleighDampingCoefficientsProcess", echo_level > 0) <<  "Solving eigen values" << std::endl;
        p_compute_eigen_values_strategy->Solve();
    }

    // The eigen values vector
    const Vector& r_eigenvalue_vector = r_current_process_info[EIGENVALUE_VECTOR];

    // We need at least 2 eigenvalues
    if (r_eigenvalue_vector.size() > 1) {
        KRATOS_INFO_IF("ComputeRayleighDampingCoefficientsProcess", echo_level > 0) <<  "The two first eigen values of the structure are:\t" << r_eigenvalue_vector[0] << " Hz\t and " << r_eigenvalue_vector[1] << " Hz" << std::endl;

        // The main damping ration
        array_1d<double, 2> damping_ratios;
        damping_ratios[0] = mParameters["damping_ratio_0"].GetDouble();
        damping_ratios[1] = mParameters["damping_ratio_0"].GetDouble();

        // We define the auxiliar alpha and beta
        double alpha = 0.0;
        double beta = 0.0;

        // Depending of the second damping ratio (if defined or not) we compute
        if (mParameters["damping_ratio_1"].GetDouble() > 0.0) {
            damping_ratios[1] = mParameters["damping_ratio_1"].GetDouble();
            BoundedMatrix<double, 2, 2> frequencies_matrix;
            frequencies_matrix(0, 0) = 0.5 / r_eigenvalue_vector[0];
            frequencies_matrix(0, 1) = 0.5 * r_eigenvalue_vector[0];
            frequencies_matrix(1, 0) = 0.5 / r_eigenvalue_vector[1];
            frequencies_matrix(1, 1) = 0.5 * r_eigenvalue_vector[1];

            // We compute the inverse
            double det;
            const BoundedMatrix<double, 2, 2> inverse_frequencies_matrix = MathUtils<double>::InvertMatrix<2>(frequencies_matrix, det);

            // We can compute now
            const array_1d<double, 2> aux_solution = prod(inverse_frequencies_matrix, damping_ratios);
            alpha = aux_solution[0];
            beta = aux_solution[1];
        } else {
            beta = (2.0 * damping_ratios[0])/(r_eigenvalue_vector[0] + r_eigenvalue_vector[1]);
            alpha = r_eigenvalue_vector[0] * r_eigenvalue_vector[1] * beta;
        }

        KRATOS_INFO_IF("ComputeRayleighDampingCoefficientsProcess", echo_level > 0) <<  "With the damping ratios of:\t" << damping_ratios[0] << "\t and " << damping_ratios[1] << "\n the RAYLEIGH_ALPHA is: " << alpha << "\tand RAYLEIGH_BETA: " << beta << std::endl;

        // We save the values
        r_current_process_info.SetValue(RAYLEIGH_ALPHA, alpha);
        r_current_process_info.SetValue(RAYLEIGH_BETA, beta);
        if (mParameters["write_on_properties"].GetBool()) {
            for (auto& r_prop : mrModelPart.rProperties()) { // The definition should be consistent and done only once
                r_prop.SetValue(RAYLEIGH_ALPHA, alpha);
                r_prop.SetValue(RAYLEIGH_BETA, beta);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters ComputeRayleighDampingCoefficientsProcess::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "echo_level"          : 0,
        "write_on_properties" : false,
        "damping_ratio_0"     : 0.0,
        "damping_ratio_1"     : -1.0,
        "eigen_system_settings" : {
            "solver_type"                : "FEASTSolver",
            "print_feast_output"         : false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem"   : true,
            "lambda_min"                 : 0.0,
            "lambda_max"                 : 4.0e5,
            "number_of_eigenvalues"      : 2,
            "search_dimension"           : 15,
            "linear_solver_settings": {
                "solver_type": "skyline_lu"
            }
        }
    })" );

    return default_parameters;
}
} // namespace Kratos
