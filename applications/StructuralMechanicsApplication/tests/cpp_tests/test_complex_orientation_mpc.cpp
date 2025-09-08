// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
// #include "includes/gid_io.h"
#include "utilities/variable_utils.h"
#include "spaces/ublas_space.h"

// Linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

// The most basic scheme (static)
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

// The most builder and solver (the block builder and solver)
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// The strategies to test
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node NodeType;

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

        // The time scheme
        typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
        typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;

        // The strategies
        typedef ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
        typedef ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;

//         void ComplexOrientationMPCGiDIODebug(ModelPart& rModelPart)
//         {
//             GidIO<> gid_io("TEST_COMPLEX_MPC", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.WriteNodalResults(DISPLACEMENT, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(LAGRANGE_DISPLACEMENT, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(FORCE_RESIDUAL, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(REACTION, rModelPart.Nodes(), label, 0);
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//         }

        void CreateComplexOrientationMPCGeometry(ModelPart& rModelPart)
        {
            rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            rModelPart.AddNodalSolutionStepVariable(REACTION);
            rModelPart.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT); // Auxiliary displacement to test complex combinations of MPC
            rModelPart.AddNodalSolutionStepVariable(FORCE_RESIDUAL); // Auxiliary displacement to test complex combinations of MPC

            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            p_elem_prop->SetValue(YOUNG_MODULUS, 210e9);
            p_elem_prop->SetValue(POISSON_RATIO, 0.3);
            p_elem_prop->SetValue(DENSITY, 1.0e4);
            array_1d<double, 3> volume_acceleration;
            volume_acceleration[0] = 0.0;
            volume_acceleration[1] = -10.0;
            volume_acceleration[2] = 0.0;
            p_elem_prop->SetValue(VOLUME_ACCELERATION, volume_acceleration);
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            // First we create the nodes
            auto p_node_1 = rModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
            auto p_node_2 = rModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
            auto p_node_3 = rModelPart.CreateNewNode( 3,  1.0 ,  1.0 , 0.0);
            auto p_node_4 = rModelPart.CreateNewNode( 4,  0.0 ,  1.0 , 0.0);

            rModelPart.CreateNewElement("SmallDisplacementElement2D4N", 1, {{1, 2, 3, 4}}, p_elem_prop);

            VariableUtils().AddDofWithReaction(DISPLACEMENT_X, REACTION_X, rModelPart);
            VariableUtils().AddDofWithReaction(DISPLACEMENT_Y, REACTION_Y, rModelPart);
            VariableUtils().AddDofWithReaction(LAGRANGE_DISPLACEMENT_X, FORCE_RESIDUAL_X, rModelPart);
            VariableUtils().AddDofWithReaction(LAGRANGE_DISPLACEMENT_Y, FORCE_RESIDUAL_Y, rModelPart);

            // Fix dofs
            p_node_1->pGetDof(DISPLACEMENT_X)->FixDof();
            p_node_1->pGetDof(DISPLACEMENT_Y)->FixDof();
            p_node_2->pGetDof(LAGRANGE_DISPLACEMENT_X)->FixDof();
//             p_node_2->pGetDof(DISPLACEMENT_Y)->FixDof();

            // The master and slave dofs
            MasterSlaveConstraint::DofPointerVectorType master_dofs(2);
            master_dofs[0] = p_node_2->pGetDof(LAGRANGE_DISPLACEMENT_X);
            master_dofs[1] = p_node_2->pGetDof(LAGRANGE_DISPLACEMENT_Y);
            MasterSlaveConstraint::DofPointerVectorType slave_dofs(2);
            slave_dofs[0] = p_node_2->pGetDof(DISPLACEMENT_X);
            slave_dofs[1] = p_node_2->pGetDof(DISPLACEMENT_Y);

            // The relation matrix and constraint one
            Matrix relation_matrix(2,2);
            relation_matrix(0, 0) = - std::sqrt(2.0)/2.0;
            relation_matrix(1, 0) =   std::sqrt(2.0)/2.0;
            relation_matrix(0, 1) =   std::sqrt(2.0)/2.0;
            relation_matrix(1, 1) =   std::sqrt(2.0)/2.0;
            Vector constant_vector = ZeroVector(2);

            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, constant_vector);

            rModelPart.CloneTimeStep(0.0);
        }

        /**
        * Checks the correct work of MPC with oriented constraints
        * Test 45º degrees constraint
        */
        KRATOS_TEST_CASE_IN_SUITE(ComplexOrientationMPC1, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("Main");
            r_model_part.SetBufferSize(2);

            CreateComplexOrientationMPCGeometry(r_model_part);

            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            BuilderAndSolverType::Pointer pbuildandsolve = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );
            SolvingStrategyType::Pointer pstrategy = SolvingStrategyType::Pointer( new ResidualBasedLinearStrategyType(r_model_part, pscheme, pbuildandsolve, true));

            pstrategy->SetEchoLevel(0);
            pstrategy->Solve();

//             // DEBUG ONLY
//             ComplexOrientationMPCGiDIODebug(r_model_part);

            KRATOS_EXPECT_LE(std::abs(r_model_part.Nodes()[2].FastGetSolutionStepValue(DISPLACEMENT_X) + 5.3651e-7), 1.0e-4);
            KRATOS_EXPECT_LE(std::abs(r_model_part.Nodes()[2].FastGetSolutionStepValue(DISPLACEMENT_Y) + 5.3651e-7), 1.0e-4);
        }


    } // namespace Testing
}  // namespace Kratos.
