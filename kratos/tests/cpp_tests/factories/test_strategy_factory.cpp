//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"

// Utility includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/strategy_factory.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        // Spaces
        using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
        using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
        using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

        /// The definition of the strategy
        using StrategyType = SolvingStrategy<SparseSpaceType,LocalSpaceType, LinearSolverType>;

        /// The definition of the factory
        using StrategyFactoryType = StrategyFactory<SparseSpaceType,LocalSpaceType, LinearSolverType>;

        /**
         * Checks if the ResidualBasedLinearStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedLinearStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "ResidualBasedLinearStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "ResidualBasedLinearStrategy");
        }

        /**
         * Checks if the ResidualBasedNewtonRaphsonStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedNewtonRaphsonStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "ResidualBasedNewtonRaphsonStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "ResidualBasedNewtonRaphsonStrategy");
        }

        /**
         * Checks if the AdaptiveResidualBasedNewtonRaphsonStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(AdaptiveResidualBasedNewtonRaphsonStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "AdaptiveResidualBasedNewtonRaphsonStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "AdaptiveResidualBasedNewtonRaphsonStrategy");
        }

        /**
         * Checks if the LineSearchStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(LineSearchStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "LineSearchStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "LineSearchStrategy");
        }

        /**
         * Checks if the ExplicitStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ExplicitStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "ExplicitStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "ExplicitStrategy");
        }

    } // namespace Testing
}  // namespace Kratos.

