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

        /// The definition of the p_builder_and_solver
        using BuilderAndSolverType = BuilderAndSolver<SparseSpaceType,LocalSpaceType, LinearSolverType>;

        /// The definition of the factory
        using BuilderAndSolverFactoryType = Factory<BuilderAndSolverType>;

        /// The definition of the explicit builder
        using ExplicitBuilderType = ExplicitBuilder<SparseSpaceType,LocalSpaceType>;

        /// The definition of the factory
        using ExplicitBuilderFactoryType = Factory<ExplicitBuilderType>;

        /// The definition of the strategy
        using StrategyType = SolvingStrategy<SparseSpaceType,LocalSpaceType, LinearSolverType>;

        /// The definition of the factory
        using StrategyFactoryType = Factory<StrategyType>;

        /// The definition of the strategy
        using ExplicitSolvingStrategyType = ExplicitSolvingStrategy<SparseSpaceType,LocalSpaceType>;

        /// The definition of the factory
        using ExplicitStrategyFactoryType = Factory<ExplicitSolvingStrategyType>;

        /// The definition of the scheme
        using SchemeType = Scheme<SparseSpaceType,LocalSpaceType>;

        /// The definition of the factory
        using SchemeFactoryType = Factory<SchemeType>;

        /// The definition of the convergence criteria
        using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType,LocalSpaceType>;

        /// The definition of the custom class
        using ConvergenceCriteriaFactoryType = Factory<ConvergenceCriteriaType>;

        /**
         * Checks if the ResidualBasedEliminationBuilderAndSolver performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedEliminationBuilderAndSolverFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"name" : "elimination_builder_and_solver"})");
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverFactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_builder_and_solver->Info(), "ResidualBasedEliminationBuilderAndSolver");
        }

        /**
         * Checks if the ResidualBasedEliminationBuilderAndSolverWithConstraints performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedEliminationBuilderAndSolverWithConstraintsFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"name" : "elimination_builder_and_solver_with_constraints"})");
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverFactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_builder_and_solver->Info(), "ResidualBasedEliminationBuilderAndSolverWithConstraints");
        }

        /**
         * Checks if the ResidualBasedBlockBuilderAndSolver performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBlockBuilderAndSolverFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"name" : "block_builder_and_solver"})");
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverFactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_builder_and_solver->Info(), "ResidualBasedBlockBuilderAndSolver");
        }

        /**
         * Checks if the ExplicitBuilder performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ExplicitBuilderFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "explicit_builder"})");
            ExplicitBuilderType::Pointer p_explicit_builder = ExplicitBuilderFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_explicit_builder->Info(), "ExplicitBuilder");
        }

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

        /**
         * Checks if the ExplicitSolvingStrategyRungeKutta4 performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ExplicitSolvingStrategyRungeKutta4Factory, KratosCoreFastSuite)
        {
            Model this_model;
            auto& r_model_part = this_model.CreateModelPart("Main");
            Parameters this_parameters = Parameters(R"({"name" : "explicit_solving_strategy_runge_kutta_4",
                                                        "explicit_builder_settings" : {
                                                            "name" : "explicit_builder"
                                                        }
                                                        })");
            ExplicitSolvingStrategyType::Pointer p_strategy = ExplicitStrategyFactoryType().Create(r_model_part, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_strategy->Info(), "ExplicitSolvingStrategyRungeKutta4");
        }

        /**
         * Checks if the ResidualBasedIncrementalUpdateStaticScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedIncrementalUpdateStaticSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "static_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedIncrementalUpdateStaticScheme");
        }

        /**
         * Checks if the ResidualBasedIncrementalUpdateStaticSchemeSlip performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedIncrementalUpdateStaticSchemeSlipFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "static_slip_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedIncrementalUpdateStaticSchemeSlip");
        }

        /**
         * Checks if the ResidualBasedBossakDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBossakDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "bossak_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedBossakDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedNewmarkDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedNewmarkDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "newmark_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedNewmarkDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedPseudoStaticDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedPseudoStaticDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "pseudo_static_scheme", "rayleigh_beta_variable" : "PRESSURE"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedPseudoStaticDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedBDFDisplacementScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBDFDisplacementSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "bdf_displacement_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedBDFDisplacementScheme");
        }

        /**
         * Checks if the ResidualBasedBDFCustomScheme performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBDFCustomSchemeFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "bdf_scheme"})");
            SchemeType::Pointer p_scheme = SchemeFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_scheme->Info(), "ResidualBasedBDFCustomScheme");
        }

        /**
         * Checks if the DisplacementCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(DisplacementCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "displacement_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "DisplacementCriteria");
        }

        /**
         * Checks if the ResidualCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "residual_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "ResidualCriteria");
        }

        /**
         * Checks if the And_Criteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(And_CriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "and_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "And_Criteria");
        }

        /**
         * Checks if the Or_Criteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(Or_CriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "or_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "Or_Criteria");
        }

        /**
         * Checks if the MixedGenericCriteria performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedGenericCriteriaFactory, KratosCoreFastSuite)
        {
            Parameters this_parameters = Parameters(R"({"name" : "mixed_generic_criteria"})");
            ConvergenceCriteriaType::Pointer p_conv_criteria = ConvergenceCriteriaFactoryType().Create(this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_conv_criteria->Info(), "MixedGenericCriteria");
        }

    } // namespace Testing
}  // namespace Kratos.

