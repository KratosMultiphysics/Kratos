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
#include "factories/register_factories.h"
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
        using StrategyFactoryType = Factory<StrategyType, ModelPart>;

        /**
         * Checks if the ResidualBasedLinearStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedLinearStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "linear_strategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "amgcl"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "static_scheme"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "elimination_builder_and_solver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_strategy->Info(), "ResidualBasedLinearStrategy");
        }

        /**
         * Checks if the ResidualBasedNewtonRaphsonStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedNewtonRaphsonStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "newton_raphson_strategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "amgcl"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "static_scheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "displacement_criteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "elimination_builder_and_solver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_strategy->Info(), "ResidualBasedNewtonRaphsonStrategy");
        }

        /**
         * Checks if the AdaptiveResidualBasedNewtonRaphsonStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(AdaptiveResidualBasedNewtonRaphsonStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "adaptative_newton_raphson_strategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "amgcl"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "static_scheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "displacement_criteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "elimination_builder_and_solver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_strategy->Info(), "AdaptiveResidualBasedNewtonRaphsonStrategy");
        }

        /**
         * Checks if the LineSearchStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(LineSearchStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "line_search_strategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "amgcl"
                                                        },
                                                        "scheme_settings" : {
                                                            "name" : "static_scheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "name" : "displacement_criteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "name" : "elimination_builder_and_solver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = StrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_strategy->Info(), "LineSearchStrategy");
        }

    } // namespace Testing
}  // namespace Kratos.

