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

            std::vector<NodeType::Pointer> indexes1({pnode1, pnode2});
            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes1});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 1, pgeom1, p_prop));
            std::vector<NodeType::Pointer> indexes2({pnode2, pnode3 });
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes2});
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

            std::vector<NodeType::Pointer> indexes1({pnode11, pnode10});
            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes1});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 1, pgeom1, p_prop));
            std::vector<NodeType::Pointer> indexes2({pnode10, pnode8 });
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes2});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 2, pgeom2, p_prop));
            std::vector<NodeType::Pointer> indexes3({pnode8, pnode6});
            GeometryType::Pointer pgeom3 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes3});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 3, pgeom3, p_prop));
            std::vector<NodeType::Pointer> indexes4({pnode6, pnode5});
            GeometryType::Pointer pgeom4 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes4});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 4, pgeom4, p_prop));
            std::vector<NodeType::Pointer> indexes5({pnode5, pnode4});
            GeometryType::Pointer pgeom5 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes5});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 5, pgeom5, p_prop));
            std::vector<NodeType::Pointer> indexes6({pnode4, pnode1});
            GeometryType::Pointer pgeom6 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes6});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 6, pgeom6, p_prop));
            std::vector<NodeType::Pointer> indexes7({pnode1, pnode2});
            GeometryType::Pointer pgeom7 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes7});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 7, pgeom7, p_prop));
            std::vector<NodeType::Pointer> indexes8({pnode2, pnode3});
            GeometryType::Pointer pgeom8 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes8});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 8, pgeom8, p_prop));
            std::vector<NodeType::Pointer> indexes9({pnode3, pnode7});
            GeometryType::Pointer pgeom9 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes9});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 9, pgeom9, p_prop));
            std::vector<NodeType::Pointer> indexes10({pnode7, pnode9});
            GeometryType::Pointer pgeom10 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes10});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 10, pgeom10, p_prop));
            std::vector<NodeType::Pointer> indexes11({pnode9, pnode11});
            GeometryType::Pointer pgeom11 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes11});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 11, pgeom11, p_prop));
            std::vector<NodeType::Pointer> indexes12({pnode10, pnode9});
            GeometryType::Pointer pgeom12 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes12});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 12, pgeom12, p_prop));
            std::vector<NodeType::Pointer> indexes13({pnode9, pnode8});
            GeometryType::Pointer pgeom13 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes13});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 13, pgeom13, p_prop));
            std::vector<NodeType::Pointer> indexes14({pnode8, pnode7});
            GeometryType::Pointer pgeom14 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes14});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 14, pgeom14, p_prop));
            std::vector<NodeType::Pointer> indexes15({pnode7, pnode6});
            GeometryType::Pointer pgeom15 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes15});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 15, pgeom15, p_prop));
            std::vector<NodeType::Pointer> indexes16({pnode6, pnode3});
            GeometryType::Pointer pgeom16 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes16});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 16, pgeom16, p_prop));
            std::vector<NodeType::Pointer> indexes17({pnode3, pnode5});
            GeometryType::Pointer pgeom17 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes17});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 17, pgeom17, p_prop));
            std::vector<NodeType::Pointer> indexes18({pnode5, pnode2});
            GeometryType::Pointer pgeom18 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes18});
            rModelPart.AddElement(Kratos::make_shared<TestBarElement>( 18, pgeom18, p_prop));
            std::vector<NodeType::Pointer> indexes19({pnode2, pnode4});
            GeometryType::Pointer pgeom19 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{indexes19});
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

//         static void DebugLHS(const SparseSpaceType::MatrixType& rA)
//         {
//             const double numeric_limit = std::numeric_limits<double>::epsilon();
//             for (int i = 0; i < rA.size1(); ++i) {
//                 for (int j = 0; j < rA.size2(); ++j) {
//                     if (std::abs(rA(i, j)) > 1.0e4 * numeric_limit) {
//                         std::cout << "            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(" << i << "," << j << ") - " << rA(i, j) << ")/rA(" << i << "," << j << ")), tolerance);" << std::endl;
//                     }
//                 }
//             }
//         }
     
        /**
         * Checks if the block builder and solver performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 6);
            KRATOS_CHECK(rA.size2() == 6);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2.069e+09)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 4.138e+09)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -2.069e+09)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -2.069e+09)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2.069e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1)/rA(5,5)), tolerance);
        }

        /**
         * Checks if the block builder and solver with constraints performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 6);
            KRATOS_CHECK(rA.size2() == 6);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2.069e+09)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 2.069e+09)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2.069e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1)/rA(5,5)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 2);
            KRATOS_CHECK(rA.size2() == 2);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 4.138e+09)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,1) - -2.069e+09)/rA(0,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,0) - -2.069e+09)/rA(1,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 2.069e+09)/rA(1,1)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver (with constraints) performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 1);
            KRATOS_CHECK(rA.size2() == 1);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2.069e+09)/rA(0,0)), tolerance);
        }

        /**
         * Checks if the block builder and solver performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 22);
            KRATOS_CHECK(rA.size2() == 22);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 7.40228e+08)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 5.98857e+08)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - 3.70114e+08)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - -1.85057e+08)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - 3.70114e+08)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1.57298e+09)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,3) - -5.55171e+08)/rA(2,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -7.40228e+08)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,5) - 3.70114e+08)/rA(2,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,9) - 3.16724e-08)/rA(2,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - -1.85057e+08)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,2) - -5.55171e+08)/rA(3,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.25748e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - 3.70114e+08)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,5) - -1.85057e+08)/rA(3,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,8) - 3.16724e-08)/rA(3,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,9) - -5.1725e+08)/rA(3,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -7.40228e+08)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - 3.70114e+08)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1.65702e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - -4.7538e+08)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,8) - -1.76565e+08)/rA(4,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,9) - -2.64848e+08)/rA(4,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,11) - 4.22299e-08)/rA(4,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,12) - -7.40228e+08)/rA(4,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,13) - 3.70114e+08)/rA(4,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,2) - 3.70114e+08)/rA(5,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,3) - -1.85057e+08)/rA(5,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - -4.7538e+08)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1.45705e+09)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,8) - -2.64848e+08)/rA(5,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,9) - -3.97272e+08)/rA(5,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,10) - 4.22299e-08)/rA(5,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,11) - -6.89667e+08)/rA(5,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,12) - 3.70114e+08)/rA(5,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,13) - -1.85057e+08)/rA(5,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 1.12703e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 7.83914e+08)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,3) - 3.16724e-08)/rA(8,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,4) - -1.76565e+08)/rA(8,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,5) - -2.64848e+08)/rA(8,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 2.24557e+09)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - 2.64848e+08)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,10) - -1.0345e+09)/rA(8,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,2) - 3.16724e-08)/rA(9,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,3) - -5.1725e+08)/rA(9,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,4) - -2.64848e+08)/rA(9,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,5) - -3.97272e+08)/rA(9,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - 2.64848e+08)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 9.14522e+08)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,5) - 4.22299e-08)/rA(10,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,8) - -1.0345e+09)/rA(10,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 2.43475e+09)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,11) - 3.65751e+08)/rA(10,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -3.65751e+08)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - -3.65751e+08)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1.0345e+09)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,4) - 4.22299e-08)/rA(11,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,5) - -6.89667e+08)/rA(11,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,10) - 3.65751e+08)/rA(11,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1.05542e+09)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - -3.65751e+08)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -3.65751e+08)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,4) - -7.40228e+08)/rA(12,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,5) - 3.70114e+08)/rA(12,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -3.65751e+08)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - -3.65751e+08)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 1.84621e+09)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -3.74477e+08)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,15) - 6.33449e-08)/rA(12,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,16) - -7.40228e+08)/rA(12,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,17) - 3.70114e+08)/rA(12,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,4) - 3.70114e+08)/rA(13,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,5) - -1.85057e+08)/rA(13,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - -3.65751e+08)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -3.65751e+08)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -3.74477e+08)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 1.77036e+09)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,14) - 6.33449e-08)/rA(13,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,15) - -1.0345e+09)/rA(13,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,16) - 3.70114e+08)/rA(13,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -1.85057e+08)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1.0345e+09)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,13) - 6.33449e-08)/rA(14,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2.80923e+09)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,15) - 3.70114e+08)/rA(14,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -7.40228e+08)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,17) - -3.70114e+08)/rA(14,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,18) - -1.0345e+09)/rA(14,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,12) - 6.33449e-08)/rA(15,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,13) - -1.0345e+09)/rA(15,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,14) - 3.70114e+08)/rA(15,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 1.21956e+09)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,16) - -3.70114e+08)/rA(15,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,17) - -1.85057e+08)/rA(15,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,12) - -7.40228e+08)/rA(16,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,13) - 3.70114e+08)/rA(16,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -7.40228e+08)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,15) - -3.70114e+08)/rA(16,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 2.22068e+09)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,17) - -3.70114e+08)/rA(16,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,19) - 1.2669e-07)/rA(16,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,20) - -7.40228e+08)/rA(16,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,21) - 3.70114e+08)/rA(16,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,12) - 3.70114e+08)/rA(17,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -1.85057e+08)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,14) - -3.70114e+08)/rA(17,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,15) - -1.85057e+08)/rA(17,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,16) - -3.70114e+08)/rA(17,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 2.62417e+09)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,18) - 1.2669e-07)/rA(17,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,19) - -2.069e+09)/rA(17,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,20) - 3.70114e+08)/rA(17,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,21) - -1.85057e+08)/rA(17,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,14) - -1.0345e+09)/rA(18,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,17) - 1.2669e-07)/rA(18,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 2.069e+09)/rA(18,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,19) - -1.2669e-07)/rA(18,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,20) - -1.0345e+09)/rA(18,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,16) - 1.2669e-07)/rA(19,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,17) - -2.069e+09)/rA(19,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,18) - -1.2669e-07)/rA(19,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,19) - 2.069e+09)/rA(19,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,16) - -7.40228e+08)/rA(20,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,17) - 3.70114e+08)/rA(20,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,18) - -1.0345e+09)/rA(20,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,20) - 1.77473e+09)/rA(20,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,21) - -3.70114e+08)/rA(20,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,16) - 3.70114e+08)/rA(21,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,17) - -1.85057e+08)/rA(21,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,20) - -3.70114e+08)/rA(21,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,21) - 1.85057e+08)/rA(21,21)), tolerance);
        }

        /**
         * Checks if the block builder and solver with constraints performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 22);
            KRATOS_CHECK(rA.size2() == 22);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 7.40228e+08)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.48622e+09)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - -1.85057e+08)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,4) - 3.70114e+08)/rA(1,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,5) - -1.85057e+08)/rA(1,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,8) - 3.16724e-08)/rA(1,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,9) - -5.1725e+08)/rA(1,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - -1.85057e+08)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1.57298e+09)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -7.40228e+08)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,5) - 3.70114e+08)/rA(2,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,9) - 3.16724e-08)/rA(2,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.25748e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,1) - 3.70114e+08)/rA(4,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -7.40228e+08)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1.65702e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - -4.7538e+08)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,8) - -1.76565e+08)/rA(4,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,9) - -2.64848e+08)/rA(4,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,11) - 4.22299e-08)/rA(4,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,12) - -7.40228e+08)/rA(4,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,13) - 3.70114e+08)/rA(4,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,1) - -1.85057e+08)/rA(5,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,2) - 3.70114e+08)/rA(5,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - -4.7538e+08)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1.45705e+09)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,8) - -2.64848e+08)/rA(5,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,9) - -3.97272e+08)/rA(5,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,10) - 4.22299e-08)/rA(5,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,11) - -6.89667e+08)/rA(5,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,12) - 3.70114e+08)/rA(5,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,13) - -1.85057e+08)/rA(5,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 1.12703e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 7.83914e+08)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,1) - 3.16724e-08)/rA(8,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,4) - -1.76565e+08)/rA(8,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,5) - -2.64848e+08)/rA(8,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 2.24557e+09)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - 2.64848e+08)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,10) - -1.0345e+09)/rA(8,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,1) - -5.1725e+08)/rA(9,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,2) - 3.16724e-08)/rA(9,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,4) - -2.64848e+08)/rA(9,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,5) - -3.97272e+08)/rA(9,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - 2.64848e+08)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 9.14522e+08)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,5) - 4.22299e-08)/rA(10,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,8) - -1.0345e+09)/rA(10,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 2.43475e+09)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,11) - 3.65751e+08)/rA(10,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -3.65751e+08)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - -3.65751e+08)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1.0345e+09)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,4) - 4.22299e-08)/rA(11,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,5) - -6.89667e+08)/rA(11,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,10) - 3.65751e+08)/rA(11,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1.05542e+09)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - -3.65751e+08)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -3.65751e+08)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,4) - -7.40228e+08)/rA(12,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,5) - 3.70114e+08)/rA(12,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -3.65751e+08)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - -3.65751e+08)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 1.84621e+09)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -3.74477e+08)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,15) - 6.33449e-08)/rA(12,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,16) - -7.40228e+08)/rA(12,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,17) - 3.70114e+08)/rA(12,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,4) - 3.70114e+08)/rA(13,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,5) - -1.85057e+08)/rA(13,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - -3.65751e+08)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -3.65751e+08)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -3.74477e+08)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 1.77036e+09)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,14) - 6.33449e-08)/rA(13,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,15) - -1.0345e+09)/rA(13,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,16) - 3.70114e+08)/rA(13,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -1.85057e+08)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1.0345e+09)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,13) - 6.33449e-08)/rA(14,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2.80923e+09)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,15) - 3.70114e+08)/rA(14,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -7.40228e+08)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,17) - -3.70114e+08)/rA(14,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,18) - -1.0345e+09)/rA(14,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,12) - 6.33449e-08)/rA(15,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,13) - -1.0345e+09)/rA(15,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,14) - 3.70114e+08)/rA(15,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 1.21956e+09)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,16) - -3.70114e+08)/rA(15,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,17) - -1.85057e+08)/rA(15,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,12) - -7.40228e+08)/rA(16,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,13) - 3.70114e+08)/rA(16,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -7.40228e+08)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,15) - -3.70114e+08)/rA(16,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 2.22068e+09)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,17) - -3.70114e+08)/rA(16,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,19) - 1.2669e-07)/rA(16,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,20) - -7.40228e+08)/rA(16,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,21) - 3.70114e+08)/rA(16,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,12) - 3.70114e+08)/rA(17,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -1.85057e+08)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,14) - -3.70114e+08)/rA(17,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,15) - -1.85057e+08)/rA(17,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,16) - -3.70114e+08)/rA(17,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 2.62417e+09)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,18) - 1.2669e-07)/rA(17,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,19) - -2.069e+09)/rA(17,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,20) - 3.70114e+08)/rA(17,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,21) - -1.85057e+08)/rA(17,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,14) - -1.0345e+09)/rA(18,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,17) - 1.2669e-07)/rA(18,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 2.069e+09)/rA(18,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,19) - -1.2669e-07)/rA(18,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,20) - -1.0345e+09)/rA(18,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,16) - 1.2669e-07)/rA(19,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,17) - -2.069e+09)/rA(19,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,18) - -1.2669e-07)/rA(19,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,19) - 2.069e+09)/rA(19,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,16) - -7.40228e+08)/rA(20,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,17) - 3.70114e+08)/rA(20,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,18) - -1.0345e+09)/rA(20,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,20) - 1.77473e+09)/rA(20,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,21) - -3.70114e+08)/rA(20,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,16) - 3.70114e+08)/rA(21,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,17) - -1.85057e+08)/rA(21,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,20) - -3.70114e+08)/rA(21,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,21) - 1.85057e+08)/rA(21,21)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 19);
            KRATOS_CHECK(rA.size2() == 19);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 5.98857e+08)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,1) - 3.70114e+08)/rA(0,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,2) - -1.85057e+08)/rA(0,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,0) - 3.70114e+08)/rA(1,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.57298e+09)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - -5.55171e+08)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - -7.40228e+08)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,4) - 3.70114e+08)/rA(1,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,0) - -1.85057e+08)/rA(2,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - -5.55171e+08)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1.25748e+09)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,3) - 3.70114e+08)/rA(2,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -1.85057e+08)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,6) - -5.1725e+08)/rA(2,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - -7.40228e+08)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,2) - 3.70114e+08)/rA(3,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.65702e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - -4.7538e+08)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,5) - -1.76565e+08)/rA(3,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,6) - -2.64848e+08)/rA(3,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,9) - -7.40228e+08)/rA(3,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,10) - 3.70114e+08)/rA(3,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,1) - 3.70114e+08)/rA(4,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -1.85057e+08)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - -4.7538e+08)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1.45705e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - -2.64848e+08)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,6) - -3.97272e+08)/rA(4,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,8) - -6.89667e+08)/rA(4,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,9) - 3.70114e+08)/rA(4,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,10) - -1.85057e+08)/rA(4,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,3) - -1.76565e+08)/rA(5,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - -2.64848e+08)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 2.24557e+09)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,6) - 2.64848e+08)/rA(5,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,7) - -1.0345e+09)/rA(5,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,2) - -5.1725e+08)/rA(6,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,3) - -2.64848e+08)/rA(6,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,4) - -3.97272e+08)/rA(6,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,5) - 2.64848e+08)/rA(6,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 9.14522e+08)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,5) - -1.0345e+09)/rA(7,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 2.43475e+09)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,8) - 3.65751e+08)/rA(7,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,9) - -3.65751e+08)/rA(7,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,10) - -3.65751e+08)/rA(7,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,11) - -1.0345e+09)/rA(7,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,4) - -6.89667e+08)/rA(8,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,7) - 3.65751e+08)/rA(8,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1.05542e+09)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - -3.65751e+08)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,10) - -3.65751e+08)/rA(8,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,3) - -7.40228e+08)/rA(9,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,4) - 3.70114e+08)/rA(9,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,7) - -3.65751e+08)/rA(9,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - -3.65751e+08)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 1.84621e+09)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,10) - -3.74477e+08)/rA(9,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,13) - -7.40228e+08)/rA(9,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,14) - 3.70114e+08)/rA(9,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,3) - 3.70114e+08)/rA(10,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,4) - -1.85057e+08)/rA(10,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,7) - -3.65751e+08)/rA(10,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,8) - -3.65751e+08)/rA(10,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,9) - -3.74477e+08)/rA(10,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 1.77036e+09)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -1.0345e+09)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - 3.70114e+08)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1.85057e+08)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,7) - -1.0345e+09)/rA(11,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 2.80923e+09)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - 3.70114e+08)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -7.40228e+08)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,14) - -3.70114e+08)/rA(11,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,15) - -1.0345e+09)/rA(11,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -1.0345e+09)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - 3.70114e+08)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 1.21956e+09)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -3.70114e+08)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,14) - -1.85057e+08)/rA(12,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,9) - -7.40228e+08)/rA(13,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - 3.70114e+08)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -7.40228e+08)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -3.70114e+08)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 2.22068e+09)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,14) - -3.70114e+08)/rA(13,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -7.40228e+08)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,18) - 3.70114e+08)/rA(13,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,9) - 3.70114e+08)/rA(14,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1.85057e+08)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,11) - -3.70114e+08)/rA(14,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,12) - -1.85057e+08)/rA(14,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,13) - -3.70114e+08)/rA(14,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2.62417e+09)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -2.069e+09)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,17) - 3.70114e+08)/rA(14,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,18) - -1.85057e+08)/rA(14,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,11) - -1.0345e+09)/rA(15,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 2.069e+09)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,17) - -1.0345e+09)/rA(15,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -2.069e+09)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 2.069e+09)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -7.40228e+08)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,14) - 3.70114e+08)/rA(17,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,15) - -1.0345e+09)/rA(17,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 1.77473e+09)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,18) - -3.70114e+08)/rA(17,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,13) - 3.70114e+08)/rA(18,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,14) - -1.85057e+08)/rA(18,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,17) - -3.70114e+08)/rA(18,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 1.85057e+08)/rA(18,18)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver (with constraints) performs correctly the resolution of the system
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
            constexpr double tolerance = 1e-4;
            KRATOS_CHECK(rA.size1() == 18);
            KRATOS_CHECK(rA.size2() == 18);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 1.48622e+09)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,1) - -1.85057e+08)/rA(0,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,2) - 3.70114e+08)/rA(0,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,3) - -1.85057e+08)/rA(0,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,5) - -5.1725e+08)/rA(0,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,0) - -1.85057e+08)/rA(1,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.57298e+09)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,2) - -7.40228e+08)/rA(1,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - 3.70114e+08)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,0) - 3.70114e+08)/rA(2,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,1) - -7.40228e+08)/rA(2,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1.65702e+09)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,3) - -4.7538e+08)/rA(2,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -1.76565e+08)/rA(2,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,5) - -2.64848e+08)/rA(2,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,8) - -7.40228e+08)/rA(2,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,9) - 3.70114e+08)/rA(2,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,0) - -1.85057e+08)/rA(3,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - 3.70114e+08)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,2) - -4.7538e+08)/rA(3,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.45705e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - -2.64848e+08)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,5) - -3.97272e+08)/rA(3,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,7) - -6.89667e+08)/rA(3,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,8) - 3.70114e+08)/rA(3,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,9) - -1.85057e+08)/rA(3,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -1.76565e+08)/rA(4,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - -2.64848e+08)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2.24557e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,5) - 2.64848e+08)/rA(4,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,6) - -1.0345e+09)/rA(4,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,0) - -5.1725e+08)/rA(5,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,2) - -2.64848e+08)/rA(5,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,3) - -3.97272e+08)/rA(5,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,4) - 2.64848e+08)/rA(5,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 9.14522e+08)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,4) - -1.0345e+09)/rA(6,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 2.43475e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,7) - 3.65751e+08)/rA(6,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,8) - -3.65751e+08)/rA(6,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,9) - -3.65751e+08)/rA(6,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,10) - -1.0345e+09)/rA(6,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,3) - -6.89667e+08)/rA(7,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,6) - 3.65751e+08)/rA(7,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 1.05542e+09)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,8) - -3.65751e+08)/rA(7,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,9) - -3.65751e+08)/rA(7,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,2) - -7.40228e+08)/rA(8,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,3) - 3.70114e+08)/rA(8,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,6) - -3.65751e+08)/rA(8,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,7) - -3.65751e+08)/rA(8,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1.84621e+09)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,9) - -3.74477e+08)/rA(8,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,12) - -7.40228e+08)/rA(8,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,13) - 3.70114e+08)/rA(8,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,2) - 3.70114e+08)/rA(9,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,3) - -1.85057e+08)/rA(9,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,6) - -3.65751e+08)/rA(9,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,7) - -3.65751e+08)/rA(9,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,8) - -3.74477e+08)/rA(9,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 1.77036e+09)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,11) - -1.0345e+09)/rA(9,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,12) - 3.70114e+08)/rA(9,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,13) - -1.85057e+08)/rA(9,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,6) - -1.0345e+09)/rA(10,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 2.80923e+09)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,11) - 3.70114e+08)/rA(10,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,12) - -7.40228e+08)/rA(10,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,13) - -3.70114e+08)/rA(10,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,14) - -1.0345e+09)/rA(10,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,9) - -1.0345e+09)/rA(11,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,10) - 3.70114e+08)/rA(11,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1.21956e+09)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,12) - -3.70114e+08)/rA(11,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,13) - -1.85057e+08)/rA(11,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,8) - -7.40228e+08)/rA(12,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,9) - 3.70114e+08)/rA(12,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,10) - -7.40228e+08)/rA(12,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,11) - -3.70114e+08)/rA(12,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 2.22068e+09)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - -3.70114e+08)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,16) - -7.40228e+08)/rA(12,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,17) - 3.70114e+08)/rA(12,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,8) - 3.70114e+08)/rA(13,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,9) - -1.85057e+08)/rA(13,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,10) - -3.70114e+08)/rA(13,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,11) - -1.85057e+08)/rA(13,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - -3.70114e+08)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 2.62417e+09)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,15) - -2.069e+09)/rA(13,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,16) - 3.70114e+08)/rA(13,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,17) - -1.85057e+08)/rA(13,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,10) - -1.0345e+09)/rA(14,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 2.069e+09)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,16) - -1.0345e+09)/rA(14,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,13) - -2.069e+09)/rA(15,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 2.069e+09)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,12) - -7.40228e+08)/rA(16,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,13) - 3.70114e+08)/rA(16,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,14) - -1.0345e+09)/rA(16,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 1.77473e+09)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,17) - -3.70114e+08)/rA(16,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,12) - 3.70114e+08)/rA(17,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,13) - -1.85057e+08)/rA(17,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,16) - -3.70114e+08)/rA(17,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 1.85057e+08)/rA(17,17)), tolerance);
        }
        
    } // namespace Testing
}  // namespace Kratos.

