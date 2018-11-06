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
#include <iomanip>

/* External includes */

/* Project includes */
#include "testing/testing.h"

/* Utility includes */
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

/* Element include */
#include "geometries/line_2d_2.h"
#include "tests/test_bar_element.h"

// Linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

// The most basic scheme (static)
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

// The builder and solvers
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"

namespace Kratos 
{
    namespace Testing 
    {
        /// Tests
        // TODO: Create test for the other components
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        
        // The direct solver
        typedef Reorderer<SparseSpaceType,  LocalSpaceType > ReordererType;
        typedef DirectSolver<SparseSpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
        typedef LinearSolver<SparseSpaceType,LocalSpaceType> LinearSolverType;
        typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;
        
        // The builder ans solver type
        typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsType;
        typedef ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;
        typedef ResidualBasedEliminationBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverWithConstraintsType;
        
        // The time scheme
        typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
        typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;
        
        /**
         * @brief It generates a truss structure with an expected solution
         */
        static inline void BasicTestBuilderAndSolverDisplacement(ModelPart& rModelPart, const bool WithConstraint = false)
        {
            rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
            rModelPart.AddNodalSolutionStepVariable(REACTION);
            rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

            NodeType::Pointer pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

            auto p_prop = rModelPart.pGetProperties(1, 0);
            p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
            p_prop->SetValue(NODAL_AREA, 0.01);

            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 2, pgeom2, p_prop));

            /// Add dof
            for (auto& node : rModelPart.Nodes()) {
                node.AddDof(DISPLACEMENT_X, REACTION_X);
                node.AddDof(DISPLACEMENT_Y, REACTION_Y);
                node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            }

            /// Initialize elements
            auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& elem : rModelPart.Elements()) {
                elem.Initialize();
                elem.InitializeSolutionStep(r_process_info);
            }

            // Set initial solution
            for (auto& node : rModelPart.Nodes()) {
                (node.FastGetSolutionStepValue(DISPLACEMENT)).clear();
                (node.FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
                (node.FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();
            }

            // Fix dofs
            for (auto& node : rModelPart.Nodes()) {
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
            }
            pnode1->Fix(DISPLACEMENT_X);

            if (WithConstraint) {
                rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode2, DISPLACEMENT_X, *pnode3, DISPLACEMENT_X, 1.0, 0.0);
            }
        }

        /**
         * @brief It generates a truss structure with an expected solution
         */
        static inline void ExtendedTestBuilderAndSolverDisplacement(ModelPart& rModelPart, const bool WithConstraint = false)
        {
            rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
            rModelPart.AddNodalSolutionStepVariable(REACTION);
            rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

            NodeType::Pointer pnode1 = rModelPart.CreateNewNode(1, 10.0, -5.0, 0.0);
            NodeType::Pointer pnode2 = rModelPart.CreateNewNode(2, 8.0, -4.0, 0.0);
            NodeType::Pointer pnode3 = rModelPart.CreateNewNode(3, 6.0, -3.0, 0.0);
            NodeType::Pointer pnode4 = rModelPart.CreateNewNode(4,10.0, 0.0, 0.0);
            NodeType::Pointer pnode5 = rModelPart.CreateNewNode(5, 8.0, 0.0, 0.0);
            NodeType::Pointer pnode6 = rModelPart.CreateNewNode(6, 6.0, 0.0, 0.0);
            NodeType::Pointer pnode7 = rModelPart.CreateNewNode(7, 4.0, -2.0, 0.0);
            NodeType::Pointer pnode8 = rModelPart.CreateNewNode(8, 4.0, 0.0, 0.0);
            NodeType::Pointer pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
            NodeType::Pointer pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);
            NodeType::Pointer pnode11 = rModelPart.CreateNewNode(11, 0.0, 0.0, 0.0);

            auto p_prop = rModelPart.pGetProperties(1, 0);
            p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
            p_prop->SetValue(NODAL_AREA, 0.01);

            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode11, pnode10})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode10, pnode8})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 2, pgeom2, p_prop));
            GeometryType::Pointer pgeom3 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode8, pnode6})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 3, pgeom3, p_prop));
            GeometryType::Pointer pgeom4 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode6, pnode5})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 4, pgeom4, p_prop));
            GeometryType::Pointer pgeom5 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode5, pnode4})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 5, pgeom5, p_prop));
            GeometryType::Pointer pgeom6 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode4, pnode1})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 6, pgeom6, p_prop));
            GeometryType::Pointer pgeom7 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 7, pgeom7, p_prop));
            GeometryType::Pointer pgeom8 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 8, pgeom8, p_prop));
            GeometryType::Pointer pgeom9 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode3, pnode7})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 9, pgeom9, p_prop));
            GeometryType::Pointer pgeom10 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode7, pnode9})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 10, pgeom10, p_prop));
            GeometryType::Pointer pgeom11 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode9, pnode11})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 11, pgeom11, p_prop));
            GeometryType::Pointer pgeom12 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode10, pnode9})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 12, pgeom12, p_prop));
            GeometryType::Pointer pgeom13 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode9, pnode8})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 13, pgeom13, p_prop));
            GeometryType::Pointer pgeom14 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode8, pnode7})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 14, pgeom14, p_prop));
            GeometryType::Pointer pgeom15 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode7, pnode6})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 15, pgeom15, p_prop));
            GeometryType::Pointer pgeom16 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode6, pnode3})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 16, pgeom16, p_prop));
            GeometryType::Pointer pgeom17 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode3, pnode5})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 17, pgeom17, p_prop));
            GeometryType::Pointer pgeom18 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode5, pnode2})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 18, pgeom18, p_prop));
            GeometryType::Pointer pgeom19 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode4})});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 19, pgeom19, p_prop));
            
            /// Add dof
            for (auto& node : rModelPart.Nodes()) {
                node.AddDof(DISPLACEMENT_X, REACTION_X);
                node.AddDof(DISPLACEMENT_Y, REACTION_Y);
                node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            }

            /// Initialize elements
            auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& elem : rModelPart.Elements()) {
                elem.Initialize();
                elem.InitializeSolutionStep(r_process_info);
            }

            // Set initial solution
            for (auto& node : rModelPart.Nodes()) {
                (node.FastGetSolutionStepValue(DISPLACEMENT)).clear();
                (node.FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
                (node.FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();
            }

            // Fix dofs
            for (auto& node : rModelPart.Nodes()) {
                node.Fix(DISPLACEMENT_Z);
            }
            pnode1->Fix(DISPLACEMENT_X);
            pnode4->Fix(DISPLACEMENT_X);
            pnode4->Fix(DISPLACEMENT_Y);

            if (WithConstraint) {
                rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode1, DISPLACEMENT_Y, *pnode2, DISPLACEMENT_Y, 1.0, 0.0);
            }
        }

        static SparseSpaceType::MatrixType BuildSystem(
            ModelPart& rModelPart,
            SchemeType::Pointer pScheme,
            BuilderAndSolverType::Pointer pBuilderAndSolver
            )
        {
            pBuilderAndSolver->SetDofSetIsInitializedFlag(false);
            pBuilderAndSolver->Clear();
            pScheme->Clear();

            SparseSpaceType::VectorPointerType pDx; /// The incremement in the solution
            SparseSpaceType::VectorPointerType pb; /// The RHS vector of the system of equations
            SparseSpaceType::MatrixPointerType pA; /// The LHS matrix of the system of equations

            pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);
            pBuilderAndSolver->SetUpSystem(rModelPart);
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, pA, pDx, pb, rModelPart);

            SparseSpaceType::MatrixType& rA  = *pA;
            SparseSpaceType::VectorType& rDx = *pDx;
            SparseSpaceType::VectorType& rb  = *pb;

            pBuilderAndSolver->InitializeSolutionStep(rModelPart, rA, rDx, rb);
            pScheme->InitializeSolutionStep(rModelPart, rA, rDx, rb);
            pScheme->InitializeNonLinIteration(rModelPart, rA, rDx, rb);

            pBuilderAndSolver->Build(pScheme, rModelPart, rA, rb);
            pBuilderAndSolver->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

            return rA;
        }

        static void DebugLHS(const SparseSpaceType::MatrixType& rA)
        {
            for (int i = 0; i < rA.size1(); ++i) {
                for (int j = 0; j < rA.size2(); ++j) {
                    if (std::abs(rA(i, j)) > 0.99) {
                        std::cout << "            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(" << i << "," << j << ") - ";
                        std::cout << std::fixed;
                        std::cout << std::setprecision(16);
                        std::cout << rA(i, j);
                        std::cout << ")/rA(" << i << "," << j << ")), tolerance);" << std::endl;
                    }
                }
            }
        }
     
        /**
         * Checks if the block builder and solver performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementBlockBuilderAndSolver, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            BasicTestBuilderAndSolverDisplacement(r_model_part);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 6);
            KRATOS_CHECK(rA.size2() == 6);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2069000000.000000000)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.000000000)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 4138000000.000000000)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -2069000000.000000000)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.000000000)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -2069000000.000000000)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2069000000.000000000)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1.000000000)/rA(5,5)), tolerance);

        }

        /**
         * Checks if the block builder and solver with constraints performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementBlockBuilderAndSolverWithConstraints, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            BasicTestBuilderAndSolverDisplacement(r_model_part, true);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverWithConstraintsType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 6);
            KRATOS_CHECK(rA.size2() == 6);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2069000000.0000000000000000)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.0000000000000000)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 2069000000.0000000000000000)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.0000000000000000)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2069000000.0000000000000000)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1.0000000000000000)/rA(5,5)), tolerance);


        }

        /**
         * Checks if the elimination builder and solver performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementEliminationBuilderAndSolver, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            BasicTestBuilderAndSolverDisplacement(r_model_part);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedEliminationBuilderAndSolverType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 2);
            KRATOS_CHECK(rA.size2() == 2);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 4138000000.000000000)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,1) - -2069000000.000000000)/rA(0,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,0) - -2069000000.000000000)/rA(1,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 2069000000.000000000)/rA(1,1)), tolerance);

        }

        /**
         * Checks if the elimination builder and solver (with constraints) performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementEliminationBuilderAndSolverWithConstraints, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            BasicTestBuilderAndSolverDisplacement(r_model_part, true);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedEliminationBuilderAndSolverWithConstraintsType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 1);
            KRATOS_CHECK(rA.size2() == 1);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2069000000.0000000000000000)/rA(0,0)), tolerance);
        }

        /**
         * Checks if the block builder and solver performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementBlockBuilderAndSolver, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            ExtendedTestBuilderAndSolverDisplacement(r_model_part);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 22);
            KRATOS_CHECK(rA.size2() == 22);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 740227943.2715302705764771)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 598856985.8178827762603760)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - 370113971.6357653141021729)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - -185056985.8178827166557312)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - 370113971.6357653141021729)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1572984379.4520018100738525)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,3) - -555170957.4536480903625488)/rA(2,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -740227943.2715302705764771)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,5) - 370113971.6357653141021729)/rA(2,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - -185056985.8178827166557312)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,2) - -555170957.4536479711532593)/rA(3,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1257477943.2715306282043457)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - 370113971.6357653141021729)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,5) - -185056985.8178827166557312)/rA(3,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,9) - -517250000.0000000000000000)/rA(3,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -740227943.2715302705764771)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - 370113971.6357653141021729)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1657021225.9261374473571777)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - -475379934.1969153881072998)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,8) - -176565339.3830768167972565)/rA(4,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,9) - -264848009.0746152102947235)/rA(4,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,12) - -740227943.2715302705764771)/rA(4,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,13) - 370113971.6357653141021729)/rA(4,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,2) - 370113971.6357652544975281)/rA(5,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,3) - -185056985.8178827166557312)/rA(5,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - -475379934.1969153881072998)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1457052651.9143548011779785)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,8) - -264848009.0746152102947235)/rA(5,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,9) - -397272013.6119228005409241)/rA(5,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,11) - -689666666.6666666269302368)/rA(5,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,12) - 370113971.6357653141021729)/rA(5,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,13) - -185056985.8178827166557312)/rA(5,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 1127028492.9089412689208984)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 783913971.6357650756835938)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,4) - -176565339.3830768167972565)/rA(8,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,5) - -264848009.0746152102947235)/rA(8,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 2245565339.3830766677856445)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - 264848009.0746151804924011)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,10) - -1034500000.0000000000000000)/rA(8,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,3) - -517250000.0000000000000000)/rA(9,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,4) - -264848009.0746152102947235)/rA(9,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,5) - -397272013.6119228005409241)/rA(9,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - 264848009.0746151804924011)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 914522013.6119227409362793)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,8) - -1034500000.0000000000000000)/rA(10,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 2434750982.5687417984008789)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,11) - 365750982.5687416195869446)/rA(10,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -365750982.5687417387962341)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - -365750982.5687416791915894)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1034500000.0000000000000000)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,5) - -689666666.6666666269302368)/rA(11,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,10) - 365750982.5687416195869446)/rA(11,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1055417649.2354083061218262)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - -365750982.5687416791915894)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -365750982.5687416195869446)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,4) - -740227943.2715302705764771)/rA(12,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,5) - 370113971.6357653141021729)/rA(12,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -365750982.5687417387962341)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - -365750982.5687416791915894)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 1846206869.1118023395538330)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -374476960.7027890682220459)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,16) - -740227943.2715302705764771)/rA(12,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,17) - 370113971.6357653141021729)/rA(12,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,4) - 370113971.6357652544975281)/rA(13,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,5) - -185056985.8178827166557312)/rA(13,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - -365750982.5687416791915894)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -365750982.5687416195869446)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -374476960.7027890682220459)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 1770364954.2045071125030518)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,15) - -1034500000.0000000000000000)/rA(13,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,16) - 370113971.6357653141021729)/rA(13,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -185056985.8178827166557312)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1034500000.0000000000000000)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2809227943.2715301513671875)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,15) - 370113971.6357650756835938)/rA(14,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -740227943.2715302705764771)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,17) - -370113971.6357651352882385)/rA(14,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,18) - -1034500000.0000000000000000)/rA(14,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,13) - -1034500000.0000000000000000)/rA(15,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,14) - 370113971.6357650756835938)/rA(15,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 1219556985.8178825378417969)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,16) - -370113971.6357651352882385)/rA(15,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,17) - -185056985.8178825676441193)/rA(15,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,12) - -740227943.2715302705764771)/rA(16,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,13) - 370113971.6357653141021729)/rA(16,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -740227943.2715302705764771)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,15) - -370113971.6357651352882385)/rA(16,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 2220683829.8145909309387207)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,17) - -370113971.6357656121253967)/rA(16,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,20) - -740227943.2715302705764771)/rA(16,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,21) - 370113971.6357653141021729)/rA(16,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,12) - 370113971.6357652544975281)/rA(17,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -185056985.8178827166557312)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,14) - -370113971.6357651352882385)/rA(17,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,15) - -185056985.8178825676441193)/rA(17,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,16) - -370113971.6357656121253967)/rA(17,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 2624170957.4536480903625488)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,19) - -2069000000.0000000000000000)/rA(17,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,20) - 370113971.6357653141021729)/rA(17,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,21) - -185056985.8178827166557312)/rA(17,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,14) - -1034500000.0000000000000000)/rA(18,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 2069000000.0000000000000000)/rA(18,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,20) - -1034500000.0000000000000000)/rA(18,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,17) - -2069000000.0000000000000000)/rA(19,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,19) - 2069000000.0000000000000000)/rA(19,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,16) - -740227943.2715302705764771)/rA(20,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,17) - 370113971.6357653141021729)/rA(20,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,18) - -1034500000.0000000000000000)/rA(20,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,20) - 1774727943.2715301513671875)/rA(20,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,21) - -370113971.6357653141021729)/rA(20,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,16) - 370113971.6357652544975281)/rA(21,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,17) - -185056985.8178827166557312)/rA(21,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,20) - -370113971.6357653141021729)/rA(21,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,21) - 185056985.8178827166557312)/rA(21,21)), tolerance);
        }

        /**
         * Checks if the block builder and solver with constraints performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementBlockBuilderAndSolverWithConstraints, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            ExtendedTestBuilderAndSolverDisplacement(r_model_part, true);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverWithConstraintsType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 22);
            KRATOS_CHECK(rA.size2() == 22);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 740227943.2715302705764771)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1486220957.4536478519439697)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - -185056985.8178826868534088)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,4) - 370113971.6357653141021729)/rA(1,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,5) - -185056985.8178827166557312)/rA(1,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,9) - -517250000.0000000000000000)/rA(1,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - -185056985.8178827762603760)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1572984379.4520018100738525)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -740227943.2715302705764771)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,5) - 370113971.6357653141021729)/rA(2,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1257477943.2715306282043457)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,1) - 370113971.6357653141021729)/rA(4,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -740227943.2715302705764771)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1657021225.9261374473571777)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - -475379934.1969153881072998)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,8) - -176565339.3830768167972565)/rA(4,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,9) - -264848009.0746152102947235)/rA(4,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,12) - -740227943.2715302705764771)/rA(4,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,13) - 370113971.6357653141021729)/rA(4,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,1) - -185056985.8178827166557312)/rA(5,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,2) - 370113971.6357652544975281)/rA(5,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - -475379934.1969153881072998)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1457052651.9143548011779785)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,8) - -264848009.0746152102947235)/rA(5,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,9) - -397272013.6119228005409241)/rA(5,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,11) - -689666666.6666666269302368)/rA(5,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,12) - 370113971.6357653141021729)/rA(5,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,13) - -185056985.8178827166557312)/rA(5,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 1127028492.9089412689208984)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 783913971.6357650756835938)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,4) - -176565339.3830768167972565)/rA(8,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,5) - -264848009.0746152102947235)/rA(8,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 2245565339.3830766677856445)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - 264848009.0746151804924011)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,10) - -1034500000.0000000000000000)/rA(8,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,1) - -517250000.0000000000000000)/rA(9,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,4) - -264848009.0746152102947235)/rA(9,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,5) - -397272013.6119228005409241)/rA(9,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - 264848009.0746151804924011)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 914522013.6119227409362793)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,8) - -1034500000.0000000000000000)/rA(10,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 2434750982.5687417984008789)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,11) - 365750982.5687416195869446)/rA(10,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -365750982.5687417387962341)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - -365750982.5687416791915894)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1034500000.0000000000000000)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,5) - -689666666.6666666269302368)/rA(11,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,10) - 365750982.5687416195869446)/rA(11,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1055417649.2354083061218262)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - -365750982.5687416791915894)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -365750982.5687416195869446)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,4) - -740227943.2715302705764771)/rA(12,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,5) - 370113971.6357653141021729)/rA(12,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -365750982.5687417387962341)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - -365750982.5687416791915894)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 1846206869.1118023395538330)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -374476960.7027890682220459)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,16) - -740227943.2715302705764771)/rA(12,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,17) - 370113971.6357653141021729)/rA(12,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,4) - 370113971.6357652544975281)/rA(13,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,5) - -185056985.8178827166557312)/rA(13,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - -365750982.5687416791915894)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -365750982.5687416195869446)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -374476960.7027890682220459)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 1770364954.2045071125030518)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,15) - -1034500000.0000000000000000)/rA(13,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,16) - 370113971.6357653141021729)/rA(13,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -185056985.8178827166557312)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1034500000.0000000000000000)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2809227943.2715301513671875)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,15) - 370113971.6357650756835938)/rA(14,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -740227943.2715302705764771)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,17) - -370113971.6357651352882385)/rA(14,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,18) - -1034500000.0000000000000000)/rA(14,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,13) - -1034500000.0000000000000000)/rA(15,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,14) - 370113971.6357650756835938)/rA(15,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 1219556985.8178825378417969)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,16) - -370113971.6357651352882385)/rA(15,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,17) - -185056985.8178825676441193)/rA(15,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,12) - -740227943.2715302705764771)/rA(16,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,13) - 370113971.6357653141021729)/rA(16,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -740227943.2715302705764771)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,15) - -370113971.6357651352882385)/rA(16,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 2220683829.8145909309387207)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,17) - -370113971.6357656121253967)/rA(16,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,20) - -740227943.2715302705764771)/rA(16,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,21) - 370113971.6357653141021729)/rA(16,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,12) - 370113971.6357652544975281)/rA(17,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -185056985.8178827166557312)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,14) - -370113971.6357651352882385)/rA(17,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,15) - -185056985.8178825676441193)/rA(17,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,16) - -370113971.6357656121253967)/rA(17,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 2624170957.4536480903625488)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,19) - -2069000000.0000000000000000)/rA(17,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,20) - 370113971.6357653141021729)/rA(17,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,21) - -185056985.8178827166557312)/rA(17,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,14) - -1034500000.0000000000000000)/rA(18,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 2069000000.0000000000000000)/rA(18,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,20) - -1034500000.0000000000000000)/rA(18,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,17) - -2069000000.0000000000000000)/rA(19,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,19) - 2069000000.0000000000000000)/rA(19,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,16) - -740227943.2715302705764771)/rA(20,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,17) - 370113971.6357653141021729)/rA(20,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,18) - -1034500000.0000000000000000)/rA(20,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,20) - 1774727943.2715301513671875)/rA(20,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,21) - -370113971.6357653141021729)/rA(20,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,16) - 370113971.6357652544975281)/rA(21,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,17) - -185056985.8178827166557312)/rA(21,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,20) - -370113971.6357653141021729)/rA(21,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,21) - 185056985.8178827166557312)/rA(21,21)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementEliminationBuilderAndSolver, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            ExtendedTestBuilderAndSolverDisplacement(r_model_part);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedEliminationBuilderAndSolverType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 19);
            KRATOS_CHECK(rA.size2() == 19);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 598856985.8178827762603760)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,1) - 370113971.6357653141021729)/rA(0,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,2) - -185056985.8178827166557312)/rA(0,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,0) - 370113971.6357653141021729)/rA(1,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1572984379.4520018100738525)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - -555170957.4536480903625488)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - -740227943.2715302705764771)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,4) - 370113971.6357653141021729)/rA(1,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,0) - -185056985.8178827166557312)/rA(2,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - -555170957.4536479711532593)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1257477943.2715306282043457)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,3) - 370113971.6357653141021729)/rA(2,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -185056985.8178827166557312)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,6) - -517250000.0000000000000000)/rA(2,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - -740227943.2715302705764771)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,2) - 370113971.6357653141021729)/rA(3,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1657021225.9261374473571777)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - -475379934.1969153881072998)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,5) - -176565339.3830768167972565)/rA(3,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,6) - -264848009.0746152102947235)/rA(3,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,9) - -740227943.2715302705764771)/rA(3,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,10) - 370113971.6357653141021729)/rA(3,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,1) - 370113971.6357652544975281)/rA(4,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -185056985.8178827166557312)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - -475379934.1969153881072998)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1457052651.9143548011779785)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - -264848009.0746152102947235)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,6) - -397272013.6119228005409241)/rA(4,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,8) - -689666666.6666666269302368)/rA(4,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,9) - 370113971.6357653141021729)/rA(4,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,10) - -185056985.8178827166557312)/rA(4,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,3) - -176565339.3830768167972565)/rA(5,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - -264848009.0746152102947235)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 2245565339.3830766677856445)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,6) - 264848009.0746151804924011)/rA(5,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,7) - -1034500000.0000000000000000)/rA(5,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,2) - -517250000.0000000000000000)/rA(6,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,3) - -264848009.0746152102947235)/rA(6,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,4) - -397272013.6119228005409241)/rA(6,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,5) - 264848009.0746151804924011)/rA(6,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 914522013.6119227409362793)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,5) - -1034500000.0000000000000000)/rA(7,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 2434750982.5687417984008789)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,8) - 365750982.5687416195869446)/rA(7,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,9) - -365750982.5687417387962341)/rA(7,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,10) - -365750982.5687416791915894)/rA(7,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,11) - -1034500000.0000000000000000)/rA(7,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,4) - -689666666.6666666269302368)/rA(8,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,7) - 365750982.5687416195869446)/rA(8,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1055417649.2354083061218262)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - -365750982.5687416791915894)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,10) - -365750982.5687416195869446)/rA(8,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,3) - -740227943.2715302705764771)/rA(9,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,4) - 370113971.6357653141021729)/rA(9,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,7) - -365750982.5687417387962341)/rA(9,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - -365750982.5687416791915894)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 1846206869.1118023395538330)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,10) - -374476960.7027890682220459)/rA(9,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,13) - -740227943.2715302705764771)/rA(9,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,14) - 370113971.6357653141021729)/rA(9,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,3) - 370113971.6357652544975281)/rA(10,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,4) - -185056985.8178827166557312)/rA(10,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,7) - -365750982.5687416791915894)/rA(10,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,8) - -365750982.5687416195869446)/rA(10,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,9) - -374476960.7027890682220459)/rA(10,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 1770364954.2045071125030518)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -1034500000.0000000000000000)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - 370113971.6357653141021729)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -185056985.8178827166557312)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,7) - -1034500000.0000000000000000)/rA(11,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 2809227943.2715301513671875)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - 370113971.6357650756835938)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -740227943.2715302705764771)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,14) - -370113971.6357651352882385)/rA(11,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,15) - -1034500000.0000000000000000)/rA(11,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -1034500000.0000000000000000)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - 370113971.6357650756835938)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 1219556985.8178825378417969)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -370113971.6357651352882385)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,14) - -185056985.8178825676441193)/rA(12,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,9) - -740227943.2715302705764771)/rA(13,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - 370113971.6357653141021729)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -740227943.2715302705764771)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -370113971.6357651352882385)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 2220683829.8145909309387207)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,14) - -370113971.6357656121253967)/rA(13,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -740227943.2715302705764771)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,18) - 370113971.6357653141021729)/rA(13,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,9) - 370113971.6357652544975281)/rA(14,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -185056985.8178827166557312)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,11) - -370113971.6357651352882385)/rA(14,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,12) - -185056985.8178825676441193)/rA(14,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,13) - -370113971.6357656121253967)/rA(14,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2624170957.4536480903625488)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -2069000000.0000000000000000)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,17) - 370113971.6357653141021729)/rA(14,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,18) - -185056985.8178827166557312)/rA(14,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,11) - -1034500000.0000000000000000)/rA(15,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 2069000000.0000000000000000)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,17) - -1034500000.0000000000000000)/rA(15,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -2069000000.0000000000000000)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 2069000000.0000000000000000)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -740227943.2715302705764771)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,14) - 370113971.6357653141021729)/rA(17,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,15) - -1034500000.0000000000000000)/rA(17,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 1774727943.2715301513671875)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,18) - -370113971.6357653141021729)/rA(17,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,13) - 370113971.6357652544975281)/rA(18,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,14) - -185056985.8178827166557312)/rA(18,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,17) - -370113971.6357653141021729)/rA(18,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 185056985.8178827166557312)/rA(18,18)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver (with constraints) performs correctly the assemble of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementEliminationBuilderAndSolverWithConstraints, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            ExtendedTestBuilderAndSolverDisplacement(r_model_part, true);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedEliminationBuilderAndSolverWithConstraintsType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

            // The solution check
            constexpr double tolerance = 1e-8;
            KRATOS_CHECK(rA.size1() == 18);
            KRATOS_CHECK(rA.size2() == 18);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 1486220957.4536478519439697)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,1) - -185056985.8178826570510864)/rA(0,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,2) - 370113971.6357653141021729)/rA(0,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,3) - -185056985.8178827166557312)/rA(0,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,5) - -517250000.0000000000000000)/rA(0,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,0) - -185056985.8178827762603760)/rA(1,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1572984379.4520018100738525)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - -740227943.2715302705764771)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - 370113971.6357653141021729)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,0) - 370113971.6357653141021729)/rA(2,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - -740227943.2715302705764771)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1657021225.9261374473571777)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,3) - -475379934.1969153881072998)/rA(2,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -176565339.3830768167972565)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,5) - -264848009.0746152102947235)/rA(2,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,8) - -740227943.2715302705764771)/rA(2,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,9) - 370113971.6357653141021729)/rA(2,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,0) - -185056985.8178827166557312)/rA(3,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - 370113971.6357652544975281)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,2) - -475379934.1969153881072998)/rA(3,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1457052651.9143548011779785)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - -264848009.0746152102947235)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,5) - -397272013.6119228005409241)/rA(3,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,7) - -689666666.6666666269302368)/rA(3,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,8) - 370113971.6357653141021729)/rA(3,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,9) - -185056985.8178827166557312)/rA(3,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -176565339.3830768167972565)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - -264848009.0746152102947235)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2245565339.3830766677856445)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - 264848009.0746151804924011)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,6) - -1034500000.0000000000000000)/rA(4,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,0) - -517250000.0000000000000000)/rA(5,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,2) - -264848009.0746152102947235)/rA(5,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,3) - -397272013.6119228005409241)/rA(5,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - 264848009.0746151804924011)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 914522013.6119227409362793)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,4) - -1034500000.0000000000000000)/rA(6,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 2434750982.5687417984008789)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,7) - 365750982.5687416195869446)/rA(6,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,8) - -365750982.5687417387962341)/rA(6,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,9) - -365750982.5687416791915894)/rA(6,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,10) - -1034500000.0000000000000000)/rA(6,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,3) - -689666666.6666666269302368)/rA(7,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,6) - 365750982.5687416195869446)/rA(7,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 1055417649.2354083061218262)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,8) - -365750982.5687416791915894)/rA(7,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,9) - -365750982.5687416195869446)/rA(7,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,2) - -740227943.2715302705764771)/rA(8,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,3) - 370113971.6357653141021729)/rA(8,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,6) - -365750982.5687417387962341)/rA(8,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,7) - -365750982.5687416791915894)/rA(8,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1846206869.1118023395538330)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - -374476960.7027890682220459)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,12) - -740227943.2715302705764771)/rA(8,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,13) - 370113971.6357653141021729)/rA(8,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,2) - 370113971.6357652544975281)/rA(9,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,3) - -185056985.8178827166557312)/rA(9,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,6) - -365750982.5687416791915894)/rA(9,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,7) - -365750982.5687416195869446)/rA(9,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - -374476960.7027890682220459)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 1770364954.2045071125030518)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,11) - -1034500000.0000000000000000)/rA(9,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,12) - 370113971.6357653141021729)/rA(9,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,13) - -185056985.8178827166557312)/rA(9,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,6) - -1034500000.0000000000000000)/rA(10,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 2809227943.2715301513671875)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,11) - 370113971.6357650756835938)/rA(10,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -740227943.2715302705764771)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - -370113971.6357651352882385)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1034500000.0000000000000000)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,9) - -1034500000.0000000000000000)/rA(11,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,10) - 370113971.6357650756835938)/rA(11,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1219556985.8178825378417969)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - -370113971.6357651352882385)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -185056985.8178825676441193)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,8) - -740227943.2715302705764771)/rA(12,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,9) - 370113971.6357653141021729)/rA(12,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -740227943.2715302705764771)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - -370113971.6357651352882385)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 2220683829.8145909309387207)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -370113971.6357656121253967)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,16) - -740227943.2715302705764771)/rA(12,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,17) - 370113971.6357653141021729)/rA(12,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,8) - 370113971.6357652544975281)/rA(13,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,9) - -185056985.8178827166557312)/rA(13,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - -370113971.6357651352882385)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -185056985.8178825676441193)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -370113971.6357656121253967)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 2624170957.4536480903625488)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,15) - -2069000000.0000000000000000)/rA(13,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,16) - 370113971.6357653141021729)/rA(13,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -185056985.8178827166557312)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1034500000.0000000000000000)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2069000000.0000000000000000)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -1034500000.0000000000000000)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,13) - -2069000000.0000000000000000)/rA(15,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 2069000000.0000000000000000)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,12) - -740227943.2715302705764771)/rA(16,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,13) - 370113971.6357653141021729)/rA(16,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -1034500000.0000000000000000)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 1774727943.2715301513671875)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,17) - -370113971.6357653141021729)/rA(16,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,12) - 370113971.6357652544975281)/rA(17,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -185056985.8178827166557312)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,16) - -370113971.6357653141021729)/rA(17,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 185056985.8178827166557312)/rA(17,17)), tolerance);
        }
        
    } // namespace Testing
}  // namespace Kratos.

