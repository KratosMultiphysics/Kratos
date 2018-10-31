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
#include "tests/test_element.h"

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
            rModelPart.AddNodalSolutionStepVariable(REACTION);
            rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

            NodeType::Pointer pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

            auto p_prop = rModelPart.pGetProperties(1, 0);
            p_prop->SetValue(DENSITY, 1500.0);
            p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
            const Variable<double>& cross_area = KratosComponents<Variable<double>>::Get("CROSS_AREA");
            p_prop->SetValue(cross_area, 0.01);
            const Variable<double>& pre_stress = KratosComponents<Variable<double>>::Get("TRUSS_PRESTRESS_PK2");
            p_prop->SetValue(pre_stress, 0.0);
            p_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw").Clone());

            std::vector<std::size_t> indexes1({1, 2});
            rModelPart.CreateNewElement("TrussElement3D2N", 1, indexes1, p_prop);
            std::vector<std::size_t> indexes2({2, 3 });
            rModelPart.CreateNewElement("TrussElement3D2N", 2, indexes2, p_prop);

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
                rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode1, DISPLACEMENT_X, *pnode2, DISPLACEMENT_X, 1.0, 0.0);
            }
        }

        /**
         * @brief It generates a truss structure with an expected solution
         */
        static inline void ExtendedTestBuilderAndSolverDisplacement(ModelPart& rModelPart, const bool WithConstraint = false)
        {
            rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
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
            p_prop->SetValue(DENSITY, 1500.0);
            p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
            const Variable<double>& cross_area = KratosComponents<Variable<double>>::Get("CROSS_AREA");
            p_prop->SetValue(cross_area, 0.01);
            const Variable<double>& pre_stress = KratosComponents<Variable<double>>::Get("TRUSS_PRESTRESS_PK2");
            p_prop->SetValue(pre_stress, 0.0);
            p_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw").Clone());

            std::vector<std::size_t> indexes1({11, 10});
            rModelPart.CreateNewElement("TrussElement3D2N", 1, indexes1, p_prop);
            std::vector<std::size_t> indexes2({10, 8 });
            rModelPart.CreateNewElement("TrussElement3D2N", 2, indexes2, p_prop);
            std::vector<std::size_t> indexes3({8, 6});
            rModelPart.CreateNewElement("TrussElement3D2N", 3, indexes3, p_prop);
            std::vector<std::size_t> indexes4({6, 5});
            rModelPart.CreateNewElement("TrussElement3D2N", 4, indexes4, p_prop);
            std::vector<std::size_t> indexes5({5, 4});
            rModelPart.CreateNewElement("TrussElement3D2N", 5, indexes5, p_prop);
            std::vector<std::size_t> indexes6({4, 1});
            rModelPart.CreateNewElement("TrussElement3D2N", 6, indexes6, p_prop);
            std::vector<std::size_t> indexes7({1, 2});
            rModelPart.CreateNewElement("TrussElement3D2N", 7, indexes7, p_prop);
            std::vector<std::size_t> indexes8({2, 3});
            rModelPart.CreateNewElement("TrussElement3D2N", 8, indexes8, p_prop);
            std::vector<std::size_t> indexes9({3, 7});
            rModelPart.CreateNewElement("TrussElement3D2N", 9, indexes9, p_prop);
            std::vector<std::size_t> indexes10({7, 9});
            rModelPart.CreateNewElement("TrussElement3D2N", 10, indexes10, p_prop);
            std::vector<std::size_t> indexes11({9, 11});
            rModelPart.CreateNewElement("TrussElement3D2N", 11, indexes11, p_prop);
            std::vector<std::size_t> indexes12({10, 9});
            rModelPart.CreateNewElement("TrussElement3D2N", 12, indexes12, p_prop);
            std::vector<std::size_t> indexes13({9, 8});
            rModelPart.CreateNewElement("TrussElement3D2N", 13, indexes13, p_prop);
            std::vector<std::size_t> indexes14({8, 7});
            rModelPart.CreateNewElement("TrussElement3D2N", 14, indexes14, p_prop);
            std::vector<std::size_t> indexes15({7, 6});
            rModelPart.CreateNewElement("TrussElement3D2N", 15, indexes15, p_prop);
            std::vector<std::size_t> indexes16({6, 3});
            rModelPart.CreateNewElement("TrussElement3D2N", 16, indexes16, p_prop);
            std::vector<std::size_t> indexes17({3, 5});
            rModelPart.CreateNewElement("TrussElement3D2N", 17, indexes17, p_prop);
            std::vector<std::size_t> indexes18({5, 2});
            rModelPart.CreateNewElement("TrussElement3D2N", 18, indexes18, p_prop);
            std::vector<std::size_t> indexes19({2, 4});
            rModelPart.CreateNewElement("TrussElement3D2N", 19, indexes19, p_prop);
            
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
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

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
            KRATOS_CHECK(rA.size1() == 9);
            KRATOS_CHECK(rA.size2() == 9);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2.069e+09)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 4.138e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,6) - -2.069e+09)/rA(3,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,3) - -2.069e+09)/rA(6,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 2.069e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 1)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1)/rA(8,8)), tolerance);
        }

        /**
         * Checks if the block builder and solver with constraints performs correctly the resolution of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementBlockBuilderAndSolverWithConstraints, KratosCoreFastSuite)
        {
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

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
            KRATOS_CHECK(rA.size1() == 9);
            KRATOS_CHECK(rA.size2() == 9);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2.069e+09)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 4.138e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 2.069e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 1)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1)/rA(8,8)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver performs correctly the resolution of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementEliminationBuilderAndSolver, KratosCoreFastSuite)
        {
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

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
        KRATOS_TEST_CASE_IN_SUITE(BasicDisplacementEliminationBuilderAndSolverWithConstraints, KratosCoreFastSuite2)
        {
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

            BasicTestBuilderAndSolverDisplacement(r_model_part, true);

            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedEliminationBuilderAndSolverWithConstraintsType(p_solver) );

            const SparseSpaceType::MatrixType& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

            // To create the solution of reference
//             DebugLHS(rA);

//             // The solution check
//             constexpr double tolerance = 1e-4;
//             KRATOS_CHECK(rA.size1() == 1);
//             KRATOS_CHECK(rA.size2() == 1);
        }

        /**
         * Checks if the block builder and solver performs correctly the resolution of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementBlockBuilderAndSolver, KratosCoreFastSuite)
        {
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

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
            KRATOS_CHECK(rA.size1() == 33);
            KRATOS_CHECK(rA.size2() == 33);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 7.40228e+08)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 5.98857e+08)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - 3.70114e+08)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,4) - -1.85057e+08)/rA(1,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - 3.70114e+08)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.57298e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,4) - -5.55171e+08)/rA(3,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,6) - -7.40228e+08)/rA(3,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,7) - 3.70114e+08)/rA(3,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,1) - -1.85057e+08)/rA(4,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,3) - -5.55171e+08)/rA(4,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1.25748e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,6) - 3.70114e+08)/rA(4,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,7) - -1.85057e+08)/rA(4,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,13) - -5.1725e+08)/rA(4,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,3) - -7.40228e+08)/rA(6,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,4) - 3.70114e+08)/rA(6,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 1.65702e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,7) - -4.7538e+08)/rA(6,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,12) - -1.76565e+08)/rA(6,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,13) - -2.64848e+08)/rA(6,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,18) - -7.40228e+08)/rA(6,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,19) - 3.70114e+08)/rA(6,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,3) - 3.70114e+08)/rA(7,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,4) - -1.85057e+08)/rA(7,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,6) - -4.7538e+08)/rA(7,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 1.45705e+09)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,12) - -2.64848e+08)/rA(7,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,13) - -3.97272e+08)/rA(7,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,16) - -6.89667e+08)/rA(7,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,18) - 3.70114e+08)/rA(7,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,19) - -1.85057e+08)/rA(7,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 1.12703e+09)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 7.83914e+08)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,6) - -1.76565e+08)/rA(12,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,7) - -2.64848e+08)/rA(12,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 2.24557e+09)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - 2.64848e+08)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,15) - -1.0345e+09)/rA(12,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,4) - -5.1725e+08)/rA(13,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,6) - -2.64848e+08)/rA(13,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,7) - -3.97272e+08)/rA(13,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - 2.64848e+08)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 9.14522e+08)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 1)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,12) - -1.0345e+09)/rA(15,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 2.43475e+09)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,16) - 3.65751e+08)/rA(15,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,18) - -3.65751e+08)/rA(15,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,19) - -3.65751e+08)/rA(15,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,21) - -1.0345e+09)/rA(15,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,7) - -6.89667e+08)/rA(16,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,15) - 3.65751e+08)/rA(16,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 1.05542e+09)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,18) - -3.65751e+08)/rA(16,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,19) - -3.65751e+08)/rA(16,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 1)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,6) - -7.40228e+08)/rA(18,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,7) - 3.70114e+08)/rA(18,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,15) - -3.65751e+08)/rA(18,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,16) - -3.65751e+08)/rA(18,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 1.84621e+09)/rA(18,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,19) - -3.74477e+08)/rA(18,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,24) - -7.40228e+08)/rA(18,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,25) - 3.70114e+08)/rA(18,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,6) - 3.70114e+08)/rA(19,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,7) - -1.85057e+08)/rA(19,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,15) - -3.65751e+08)/rA(19,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,16) - -3.65751e+08)/rA(19,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,18) - -3.74477e+08)/rA(19,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,19) - 1.77036e+09)/rA(19,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,22) - -1.0345e+09)/rA(19,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,24) - 3.70114e+08)/rA(19,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,25) - -1.85057e+08)/rA(19,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,20) - 1)/rA(20,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,15) - -1.0345e+09)/rA(21,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,21) - 2.80923e+09)/rA(21,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,22) - 3.70114e+08)/rA(21,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,24) - -7.40228e+08)/rA(21,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,25) - -3.70114e+08)/rA(21,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,27) - -1.0345e+09)/rA(21,27)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,19) - -1.0345e+09)/rA(22,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,21) - 3.70114e+08)/rA(22,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,22) - 1.21956e+09)/rA(22,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,24) - -3.70114e+08)/rA(22,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,25) - -1.85057e+08)/rA(22,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(23,23) - 1)/rA(23,23)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,18) - -7.40228e+08)/rA(24,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,19) - 3.70114e+08)/rA(24,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,21) - -7.40228e+08)/rA(24,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,22) - -3.70114e+08)/rA(24,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,24) - 2.22068e+09)/rA(24,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,25) - -3.70114e+08)/rA(24,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,30) - -7.40228e+08)/rA(24,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,31) - 3.70114e+08)/rA(24,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,18) - 3.70114e+08)/rA(25,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,19) - -1.85057e+08)/rA(25,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,21) - -3.70114e+08)/rA(25,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,22) - -1.85057e+08)/rA(25,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,24) - -3.70114e+08)/rA(25,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,25) - 2.62417e+09)/rA(25,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,28) - -2.069e+09)/rA(25,28)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,30) - 3.70114e+08)/rA(25,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,31) - -1.85057e+08)/rA(25,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(26,26) - 1)/rA(26,26)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(27,21) - -1.0345e+09)/rA(27,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(27,27) - 2.069e+09)/rA(27,27)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(27,30) - -1.0345e+09)/rA(27,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(28,25) - -2.069e+09)/rA(28,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(28,28) - 2.069e+09)/rA(28,28)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(29,29) - 1)/rA(29,29)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,24) - -7.40228e+08)/rA(30,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,25) - 3.70114e+08)/rA(30,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,27) - -1.0345e+09)/rA(30,27)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,30) - 1.77473e+09)/rA(30,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,31) - -3.70114e+08)/rA(30,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,24) - 3.70114e+08)/rA(31,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,25) - -1.85057e+08)/rA(31,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,30) - -3.70114e+08)/rA(31,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,31) - 1.85057e+08)/rA(31,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(32,32) - 1)/rA(32,32)), tolerance);
        }

        /**
         * Checks if the block builder and solver with constraints performs correctly the resolution of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementBlockBuilderAndSolverWithConstraints, KratosCoreFastSuite)
        {
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

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
            KRATOS_CHECK(rA.size1() == 33);
            KRATOS_CHECK(rA.size2() == 33);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 7.40228e+08)/rA(0,0)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.48622e+09)/rA(1,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,3) - -1.85057e+08)/rA(1,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,6) - 3.70114e+08)/rA(1,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,7) - -1.85057e+08)/rA(1,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,13) - -5.1725e+08)/rA(1,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 1)/rA(2,2)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,1) - -1.85057e+08)/rA(3,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.57298e+09)/rA(3,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,6) - -7.40228e+08)/rA(3,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,7) - 3.70114e+08)/rA(3,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 1.25748e+09)/rA(4,4)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1)/rA(5,5)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,1) - 3.70114e+08)/rA(6,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,3) - -7.40228e+08)/rA(6,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,6) - 1.65702e+09)/rA(6,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,7) - -4.7538e+08)/rA(6,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,12) - -1.76565e+08)/rA(6,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,13) - -2.64848e+08)/rA(6,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,18) - -7.40228e+08)/rA(6,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(6,19) - 3.70114e+08)/rA(6,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,1) - -1.85057e+08)/rA(7,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,3) - 3.70114e+08)/rA(7,3)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,6) - -4.7538e+08)/rA(7,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,7) - 1.45705e+09)/rA(7,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,12) - -2.64848e+08)/rA(7,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,13) - -3.97272e+08)/rA(7,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,16) - -6.89667e+08)/rA(7,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,18) - 3.70114e+08)/rA(7,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(7,19) - -1.85057e+08)/rA(7,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(8,8) - 1)/rA(8,8)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(9,9) - 1.12703e+09)/rA(9,9)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(10,10) - 7.83914e+08)/rA(10,10)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(11,11) - 1)/rA(11,11)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,6) - -1.76565e+08)/rA(12,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,7) - -2.64848e+08)/rA(12,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,12) - 2.24557e+09)/rA(12,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,13) - 2.64848e+08)/rA(12,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(12,15) - -1.0345e+09)/rA(12,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,1) - -5.1725e+08)/rA(13,1)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,6) - -2.64848e+08)/rA(13,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,7) - -3.97272e+08)/rA(13,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,12) - 2.64848e+08)/rA(13,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(13,13) - 9.14522e+08)/rA(13,13)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(14,14) - 1)/rA(14,14)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,12) - -1.0345e+09)/rA(15,12)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,15) - 2.43475e+09)/rA(15,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,16) - 3.65751e+08)/rA(15,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,18) - -3.65751e+08)/rA(15,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,19) - -3.65751e+08)/rA(15,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(15,21) - -1.0345e+09)/rA(15,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,7) - -6.89667e+08)/rA(16,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,15) - 3.65751e+08)/rA(16,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,16) - 1.05542e+09)/rA(16,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,18) - -3.65751e+08)/rA(16,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(16,19) - -3.65751e+08)/rA(16,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(17,17) - 1)/rA(17,17)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,6) - -7.40228e+08)/rA(18,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,7) - 3.70114e+08)/rA(18,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,15) - -3.65751e+08)/rA(18,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,16) - -3.65751e+08)/rA(18,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,18) - 1.84621e+09)/rA(18,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,19) - -3.74477e+08)/rA(18,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,24) - -7.40228e+08)/rA(18,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(18,25) - 3.70114e+08)/rA(18,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,6) - 3.70114e+08)/rA(19,6)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,7) - -1.85057e+08)/rA(19,7)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,15) - -3.65751e+08)/rA(19,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,16) - -3.65751e+08)/rA(19,16)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,18) - -3.74477e+08)/rA(19,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,19) - 1.77036e+09)/rA(19,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,22) - -1.0345e+09)/rA(19,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,24) - 3.70114e+08)/rA(19,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(19,25) - -1.85057e+08)/rA(19,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(20,20) - 1)/rA(20,20)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,15) - -1.0345e+09)/rA(21,15)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,21) - 2.80923e+09)/rA(21,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,22) - 3.70114e+08)/rA(21,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,24) - -7.40228e+08)/rA(21,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,25) - -3.70114e+08)/rA(21,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(21,27) - -1.0345e+09)/rA(21,27)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,19) - -1.0345e+09)/rA(22,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,21) - 3.70114e+08)/rA(22,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,22) - 1.21956e+09)/rA(22,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,24) - -3.70114e+08)/rA(22,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(22,25) - -1.85057e+08)/rA(22,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(23,23) - 1)/rA(23,23)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,18) - -7.40228e+08)/rA(24,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,19) - 3.70114e+08)/rA(24,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,21) - -7.40228e+08)/rA(24,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,22) - -3.70114e+08)/rA(24,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,24) - 2.22068e+09)/rA(24,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,25) - -3.70114e+08)/rA(24,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,30) - -7.40228e+08)/rA(24,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(24,31) - 3.70114e+08)/rA(24,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,18) - 3.70114e+08)/rA(25,18)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,19) - -1.85057e+08)/rA(25,19)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,21) - -3.70114e+08)/rA(25,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,22) - -1.85057e+08)/rA(25,22)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,24) - -3.70114e+08)/rA(25,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,25) - 2.62417e+09)/rA(25,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,28) - -2.069e+09)/rA(25,28)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,30) - 3.70114e+08)/rA(25,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(25,31) - -1.85057e+08)/rA(25,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(26,26) - 1)/rA(26,26)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(27,21) - -1.0345e+09)/rA(27,21)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(27,27) - 2.069e+09)/rA(27,27)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(27,30) - -1.0345e+09)/rA(27,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(28,25) - -2.069e+09)/rA(28,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(28,28) - 2.069e+09)/rA(28,28)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(29,29) - 1)/rA(29,29)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,24) - -7.40228e+08)/rA(30,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,25) - 3.70114e+08)/rA(30,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,27) - -1.0345e+09)/rA(30,27)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,30) - 1.77473e+09)/rA(30,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(30,31) - -3.70114e+08)/rA(30,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,24) - 3.70114e+08)/rA(31,24)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,25) - -1.85057e+08)/rA(31,25)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,30) - -3.70114e+08)/rA(31,30)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(31,31) - 1.85057e+08)/rA(31,31)), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((rA(32,32) - 1)/rA(32,32)), tolerance);
        }

        /**
         * Checks if the elimination builder and solver performs correctly the resolution of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(ExtendedDisplacementEliminationBuilderAndSolver, KratosCoreFastSuite)
        {
            if (!KratosComponents<Element>::Has("TrussElement3D2N")) {
                std::cout << "Please compile the StructuralMechanicsApplication in order to run this test" << std::endl;
                return void();
            }

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
        
    } // namespace Testing
}  // namespace Kratos.

