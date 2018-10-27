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
       
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        
        /// The definition of the strategy
        typedef SolvingStrategy<SparseSpaceType,LocalSpaceType, LinearSolverType> StrategyType;

        /// The definition of the factory
        typedef StrategyFactory<SparseSpaceType,LocalSpaceType, LinearSolverType> FactoryType;
     
        /**
         * Checks if the ResidualBasedLinearStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedLinearStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"strategy_type" : "ResidualBasedLinearStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "scheme_type" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "convergence_criterion" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "builder_and_solver_type" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = FactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "ResidualBasedLinearStrategy");
        }

        /**
         * Checks if the ResidualBasedNewtonRaphsonStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedNewtonRaphsonStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"strategy_type" : "ResidualBasedNewtonRaphsonStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "scheme_type" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "convergence_criterion" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "builder_and_solver_type" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = FactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "ResidualBasedNewtonRaphsonStrategy");
        }

        /**
         * Checks if the AdaptiveResidualBasedNewtonRaphsonStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(AdaptiveResidualBasedNewtonRaphsonStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"strategy_type" : "AdaptiveResidualBasedNewtonRaphsonStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "scheme_type" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "convergence_criterion" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "builder_and_solver_type" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = FactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "AdaptiveResidualBasedNewtonRaphsonStrategy");
        }

        /**
         * Checks if the LineSearchStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(LineSearchStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"strategy_type" : "LineSearchStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "scheme_type" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "convergence_criterion" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "builder_and_solver_type" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = FactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "LineSearchStrategy");
        }

        /**
         * Checks if the ExplicitStrategy performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ExplicitStrategyFactory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"strategy_type" : "ExplicitStrategy",
                                                        "linear_solver_settings" : {
                                                            "solver_type" : "AMGCL"
                                                        },
                                                        "scheme_settings" : {
                                                            "scheme_type" : "ResidualBasedIncrementalUpdateStaticScheme"
                                                        },
                                                        "convergence_criteria_settings" : {
                                                            "convergence_criterion" : "DisplacementCriteria"
                                                        },
                                                        "builder_and_solver_settings" : {
                                                            "builder_and_solver_type" : "ResidualBasedEliminationBuilderAndSolver"
                                                        }
                                                        })");
            StrategyType::Pointer p_strategy = FactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_strategy->Info()).c_str(), "ExplicitStrategy");
        }
        
    } // namespace Testing
}  // namespace Kratos.

