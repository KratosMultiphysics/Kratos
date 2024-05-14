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
//                   Alejandro Cornejo
//

// System includes
#include <limits>

/* External includes */

/* Project includes */
#include "testing/testing.h"

/* Utility includes */
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

/* Element include */
#include "tests/test_utilities/test_element.h"
#include "tests/test_utilities/test_constitutive_law.h"

// Linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

// The most basic scheme (static)
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

// The most extended criteria (residual criteria)
#include "solving_strategies/convergencecriterias/residual_criteria.h"

// The most builder and solver (the block builder and solver)
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// The strategies to test
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "solving_strategies/strategies/arc_length_strategy.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests
        // NOTE: The strategies test many things simulataneously
        // TODO: Create test for the other components
        typedef Node NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        // The direct solver
        typedef Reorderer<SparseSpaceType,  LocalSpaceType > ReordererType;
        typedef DirectSolver<SparseSpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
        typedef LinearSolver<SparseSpaceType,LocalSpaceType> LinearSolverType;
        typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

        // The convergence criteria type
        typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
        typedef ResidualCriteria< SparseSpaceType, LocalSpaceType > ResidualCriteriaType;

        // The builder ans solver type
        typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;

        // The time scheme
        typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
        typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;

        // The strategies
        typedef ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
        typedef ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;
        typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;
        typedef ArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ArcLengthStrategyType;
        typedef LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >LineSearchStrategyType;

        // Dof arrays
        typedef ModelPart::DofsArrayType DofsArrayType;
        typedef TestElement::ResidualType ResidualType;

        static inline DofsArrayType BasicTestStrategyDisplacement(
            ModelPart& ModelPart,
            const ResidualType ThisResidualType
            )
        {
            ModelPart.SetBufferSize(3);

            ModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            ModelPart.AddNodalSolutionStepVariable(REACTION);
            ModelPart.AddNodalSolutionStepVariable(VELOCITY);
            ModelPart.AddNodalSolutionStepVariable(ACCELERATION);

            Properties::Pointer p_prop = ModelPart.CreateNewProperties(0);
            TestConstitutiveLaw r_clone_cl = TestConstitutiveLaw();
            p_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            NodeType::Pointer pnode = ModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            std::vector<NodeType::Pointer> geom(1);
            geom[0] = pnode;
            GeometryType::Pointer pgeom = Kratos::make_shared<GeometryType>(PointerVector<NodeType>{geom});
            Element::Pointer pelem = Kratos::make_intrusive<TestElement>(1, pgeom, p_prop, ThisResidualType);
            ModelPart.AddElement(pelem);

            pnode->AddDof(DISPLACEMENT_X, REACTION_X);
            pnode->AddDof(DISPLACEMENT_Y, REACTION_Y);
            pnode->AddDof(DISPLACEMENT_Z, REACTION_Z);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(3);
            DoF.push_back(pnode->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode->pGetDof(DISPLACEMENT_Z));

            // Set initial solution
            (pnode->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode->FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
            (pnode->FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            Doftemp.Sort();

            return Doftemp;
        }

        /**
         * Checks if the Linear strategy performs correctly the resolution of the system
         */

        KRATOS_TEST_CASE_IN_SUITE(DisplacementLinearStrategy, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance = 1e-6;

            ModelPart& model_part = current_model.CreateModelPart("Main");

            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer pbuildandsolve = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );

            SolvingStrategyType::Pointer pstrategy = SolvingStrategyType::Pointer( new ResidualBasedLinearStrategyType(model_part, pscheme, pbuildandsolve, true));

            DofsArrayType Doftemp = BasicTestStrategyDisplacement(model_part, ResidualType::LINEAR);

            NodeType::Pointer pnode = model_part.pGetNode(1);

            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;
            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;

                model_part.CloneTimeStep(time);

                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                for (auto it= Doftemp.begin(); it!= Doftemp.end(); it++) {
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepValue() - 1.0), tolerance);
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepReactionValue()), tolerance);
                }
            }
        }

        /**
         * Checks if the Newton Rapshon strategy performs correctly the resolution of the system
         */

        KRATOS_TEST_CASE_IN_SUITE(DisplacementNRStrategy, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance = 1e-6;

            ModelPart& model_part = current_model.CreateModelPart("Main");

            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            ConvergenceCriteriaType::Pointer pcriteria = ConvergenceCriteriaType::Pointer( new ResidualCriteriaType(1.0e-4, 1.0e-9) );
            BuilderAndSolverType::Pointer pbuildandsolve = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );

            SolvingStrategyType::Pointer pstrategy = SolvingStrategyType::Pointer( new ResidualBasedNewtonRaphsonStrategyType(model_part, pscheme, pcriteria, pbuildandsolve, 10, true));

            DofsArrayType Doftemp = BasicTestStrategyDisplacement(model_part, ResidualType::NON_LINEAR);

            NodeType::Pointer pnode = model_part.pGetNode(1);

            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;
            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;

                model_part.CloneTimeStep(time);

                array_1d<double, 3> init_vector;
                init_vector[0] = 0.5;
                init_vector[1] = 0.5;
                init_vector[2] = 0.5;
                pnode->FastGetSolutionStepValue(DISPLACEMENT) = init_vector;

                pcriteria->SetEchoLevel(0);
                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                for (auto it= Doftemp.begin(); it!= Doftemp.end(); it++) {
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepValue() - 1.0), tolerance);
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepReactionValue()), tolerance);
                }
            }
        }


        /**
         * Checks if the Nonconvereged solutions of the Newton Rapshon strategy performs correctly
         */
        KRATOS_TEST_CASE_IN_SUITE(NonconvergedSolutionsNRStrategy, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance_residual_criteria = 1e-15;  // with this tolerance, I expect the strategy to reach the maximum number of iterations, i.e. 10


            ModelPart& model_part = current_model.CreateModelPart("Main");

            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            ConvergenceCriteriaType::Pointer pcriteria = ConvergenceCriteriaType::Pointer( new ResidualCriteriaType(tolerance_residual_criteria, tolerance_residual_criteria) );
            BuilderAndSolverType::Pointer pbuildandsolve = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );

            ResidualBasedNewtonRaphsonStrategyType::Pointer pstrategy = ResidualBasedNewtonRaphsonStrategyType::Pointer( new ResidualBasedNewtonRaphsonStrategyType(model_part, pscheme, pcriteria, pbuildandsolve, 10, true));

            DofsArrayType Doftemp = BasicTestStrategyDisplacement(model_part, ResidualType::NON_LINEAR);

            NodeType::Pointer pnode = model_part.pGetNode(1);

            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;

            pstrategy-> SetUpNonconvergedSolutionsGathering();

            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;

                model_part.CloneTimeStep(time);

                array_1d<double, 3> init_vector;
                init_vector[0] = 0.5;
                init_vector[1] = 0.5;
                init_vector[2] = 0.5;
                pnode->FastGetSolutionStepValue(DISPLACEMENT) = init_vector;

                pcriteria->SetEchoLevel(0);
                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                auto [NonconveregedSolutionsMatrix,Dofs]= pstrategy->GetNonconvergedSolutions();

                unsigned int expected_rows = Dofs.size();
                unsigned int expected_cols = model_part.GetProcessInfo()[NL_ITERATION_NUMBER] + 1; //+1 because zeroth is included

                KRATOS_CHECK_EQUAL(NonconveregedSolutionsMatrix.size1(), expected_rows);
                KRATOS_CHECK_EQUAL(NonconveregedSolutionsMatrix.size2(), expected_cols);


            }
        }


        KRATOS_TEST_CASE_IN_SUITE(LineSearchStrategy, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance = 1e-6;

            ModelPart& model_part = current_model.CreateModelPart("Main");

            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            ConvergenceCriteriaType::Pointer pcriteria = ConvergenceCriteriaType::Pointer( new ResidualCriteriaType(1.0e-4, 1.0e-9) );
            BuilderAndSolverType::Pointer pbuildandsolve = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );
            Parameters settings(R"({
                "max_iteration"              : 10,
                "compute_reactions"        : true,
                "max_line_search_iterations" : 5,
                "first_alpha_value"          : 0.5,
                "second_alpha_value"         : 1.0,
                "min_alpha"                  : 0.1,
                "max_alpha"                  : 2.0,
                "line_search_tolerance"      : 0.5
            })");

            SolvingStrategyType::Pointer pstrategy = SolvingStrategyType::Pointer( new LineSearchStrategyType(model_part, pscheme, pcriteria, pbuildandsolve, settings));

            DofsArrayType Doftemp = BasicTestStrategyDisplacement(model_part, ResidualType::NON_LINEAR);

            NodeType::Pointer pnode = model_part.pGetNode(1);

            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;
            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;

                model_part.CloneTimeStep(time);

                array_1d<double, 3> init_vector;
                init_vector[0] = 0.5;
                init_vector[1] = 0.5;
                init_vector[2] = 0.5;
                pnode->FastGetSolutionStepValue(DISPLACEMENT) = init_vector;

                pcriteria->SetEchoLevel(0);
                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                for (auto it= Doftemp.begin(); it!= Doftemp.end(); it++) {
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepValue() - 1.0), tolerance);
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepReactionValue()), tolerance);
                }
            }
        }

        /**
         * Checks if the Newton Rapshon strategy performs correctly the resolution of the system
         */

        KRATOS_TEST_CASE_IN_SUITE(DisplacementArcLengthStrategy, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance = 1e-6;

            ModelPart& model_part = current_model.CreateModelPart("Main");

            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            ConvergenceCriteriaType::Pointer pcriteria = ConvergenceCriteriaType::Pointer( new ResidualCriteriaType(1.0e-4, 1.0e-9) );
            BuilderAndSolverType::Pointer pbuildandsolve = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );
            Parameters settings(R"({

                "max_iteration"            : 30,
                "compute_reactions"        : true,
                "reform_dofs_at_each_step" : true,
                "move_mesh_flag"           : true,
                "desired_iterations"       : 4,
                "max_radius_factor"        : 10.0,
                "min_radius_factor"        : 0.1,
                "loads_sub_model_part_list": [],
                "loads_variable_list"      : []
            })");

            SolvingStrategyType::Pointer pstrategy = SolvingStrategyType::Pointer( new ArcLengthStrategyType(model_part, pscheme, pcriteria, pbuildandsolve, settings));

            DofsArrayType Doftemp = BasicTestStrategyDisplacement(model_part, ResidualType::NON_LINEAR);

            NodeType::Pointer pnode = model_part.pGetNode(1);

            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;
            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;

                model_part.CloneTimeStep(time);

                array_1d<double, 3> init_vector;
                init_vector[0] = 0.5;
                init_vector[1] = 0.5;
                init_vector[2] = 0.5;
                pnode->FastGetSolutionStepValue(DISPLACEMENT) = init_vector;

                pcriteria->SetEchoLevel(0);
                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                for (auto it= Doftemp.begin(); it!= Doftemp.end(); it++) {
                    KRATOS_EXPECT_LE(std::abs(it->GetSolutionStepValue() - 1.25), tolerance);
                }
            }
        }

    } // namespace Testing
}  // namespace Kratos.

