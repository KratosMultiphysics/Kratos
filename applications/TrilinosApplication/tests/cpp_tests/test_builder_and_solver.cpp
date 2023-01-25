//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>
#include <iomanip>

// External includes

/* Trilinos includes */
#include "Epetra_FEVector.h"

// Project includes
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "trilinos_space.h"
#include "containers/model.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/model_part_communicator_utilities.h"

/* Element include */
#include "geometries/line_2d_2.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_bar_element.h"

/* Linear solvers */
#include "linear_solvers/linear_solver.h"
#include "external_includes/amgcl_mpi_solver.h"

// The most basic scheme (static)
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

// The builder and solvers
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_elimination_builder_and_solver.h"

namespace Kratos::Testing
{
    /// Initial definitons
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    // Linear solvers Trilinos definitions
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISolverType;

    // Builder and solvers definition
    typedef BuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;
    typedef TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBlockBuilderAndSolverType;
    typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosResidualBasedEliminationBuilderAndSolverType;

    // The time scheme
    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType >  TrilinosSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosResidualBasedIncrementalUpdateStaticSchemeType;

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void BasicTestBuilderAndSolverDisplacement(
        ModelPart& rModelPart,
        const bool WithConstraint = false,
        const bool AdditionalNode = false
        )
    {
        rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);
        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(REACTION);
        rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Create properties
        auto p_prop = rModelPart.CreateNewProperties(1, 0);
        p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
        p_prop->SetValue(NODAL_AREA, 0.01);

        // MPI data
        const DataCommunicator& r_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int rank =  r_comm.Rank();
        const int world_size = r_comm.Size();

        // Initially everything in one partition
        if (rank == 0) {
            // Create nodes
            auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            auto pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

            // Create elements
            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
        }

        /// Add PARTITION_INDEX
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
        }

        /// Add dof
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X, REACTION_X);
            r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        }

        /// Initialize elements
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements()) {
            r_elem.Initialize(r_process_info);
            r_elem.InitializeSolutionStep(r_process_info);
        }

        // Set initial solution
        for (auto& r_node : rModelPart.Nodes()) {
            (r_node.FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (r_node.FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
            (r_node.FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();
        }

        // Fix dofs
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.Fix(DISPLACEMENT_Y);
            r_node.Fix(DISPLACEMENT_Z);
        }
        // Fix X in first node
        if (rModelPart.HasNode(1)) {
            auto pnode = rModelPart.pGetNode(1);
            pnode->Fix(DISPLACEMENT_X);
        }

        // TODO: Fix this!
        // if (WithConstraint) {
        //     rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode2, DISPLACEMENT_X, *pnode3, DISPLACEMENT_X, 1.0, 0.0);
        //     if (AdditionalNode) {
        //         rModelPart.CreateNewNode(4, 3.0, 0.0, 0.0);
        //         rModelPart.CreateNewNode(5, 4.0, 0.0, 0.0);
        //         auto pnode4 = rModelPart.pGetNode(4);
        //         rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, *pnode3, DISPLACEMENT_X, *pnode4, DISPLACEMENT_X, 1.0, 0.0);
        //         auto pnode5 = rModelPart.pGetNode(5);
        //         rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, *pnode3, DISPLACEMENT_X, *pnode5, DISPLACEMENT_X, 1.0, 0.0);
        //     }
        // }
    }

    static TrilinosSparseSpaceType::MatrixType BuildSystem(
        ModelPart& rModelPart,
        TrilinosSchemeType::Pointer pScheme,
        TrilinosBuilderAndSolverType::Pointer pBuilderAndSolver
        )
    {
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        pBuilderAndSolver->Clear();
        pScheme->Clear();

        TrilinosSparseSpaceType::VectorPointerType pDx; /// The incremement in the solution
        TrilinosSparseSpaceType::VectorPointerType pb; /// The RHS vector of the system of equations
        TrilinosSparseSpaceType::MatrixPointerType pA; /// The LHS matrix of the system of equations

        pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);
        pBuilderAndSolver->SetUpSystem(rModelPart);
        pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, pA, pDx, pb, rModelPart);

        TrilinosSparseSpaceType::MatrixType& rA  = *pA;
        TrilinosSparseSpaceType::VectorType& rDx = *pDx;
        TrilinosSparseSpaceType::VectorType& rb  = *pb;

        pBuilderAndSolver->InitializeSolutionStep(rModelPart, rA, rDx, rb);
        pScheme->InitializeSolutionStep(rModelPart, rA, rDx, rb);
        pScheme->InitializeNonLinIteration(rModelPart, rA, rDx, rb);

        pBuilderAndSolver->Build(pScheme, rModelPart, rA, rb);
        // TODO: Activate when MPC work
        // if(rModelPart.MasterSlaveConstraints().size() != 0) {
        //     pBuilderAndSolver->ApplyConstraints(pScheme, rModelPart, rA, rb);
        // }
        pBuilderAndSolver->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        return rA;
    }

    // NOTE: Adapt this to MPI
    static void DebugLHS(
        const TrilinosSparseSpaceType::MatrixType& rA,
        const DataCommunicator& rCommunicator
        )
    {
        std::cout << "KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), " << rA.NumGlobalRows() << ");";
        std::cout << "KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), " << rA.NumGlobalCols() << ");";
        std::cout << "KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), " << rA.NumGlobalNonzeros() << ");";
        // const int rank =  rCommunicator.Rank();
        // for (std::size_t i = 0; i < rA.size1(); ++i) {
        //     for (std::size_t j = 0; j < rA.size2(); ++j) {
        //         if (std::abs(rA(i, j)) > 0.99) {
        //             std::cout << "            if (rank == " << rank << ") KRATOS_CHECK_RELATIVE_NEAR(rA(" << i << "," << j << "), ";
        //             std::cout << std::fixed;
        //             std::cout << std::setprecision(16);
        //             std::cout << rA(i, j);
        //             std::cout << ", tolerance);" << std::endl;
        //         }
        //     }
        // }
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system
    */
    KRATOS_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
    {
        // The base model part
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

        // The data communicator
        const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
        const int rank =  r_comm.Rank();
        const int world_size = r_comm.Size();

        // Generate Epetra coomunicator
        KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        // Set MPI coomunicator
        ModelPartCommunicatorUtilities::SetMPICommunicator(r_model_part, r_comm);

        // Basic build
        BasicTestBuilderAndSolverDisplacement(r_model_part);

        // Compute communicaton plan and fill communicator meshes correctly
        ParallelFillCommunicator(r_model_part, r_comm).Execute();

        // Create the solvers and things required
        auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
        auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
        auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, 15, p_solver) );

        const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

        // // To create the solution of reference
        // DebugLHS(rA);

        // The solution check
        constexpr double tolerance = 1e-8;
        KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 28);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(1,1), 1.0, tolerance);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(2,2), 4138000000.0, tolerance);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(2,4), -2069000000.0, tolerance);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(3,3), 1.0, tolerance);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(4,2), -2069000000.0, tolerance);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(4,4), 2069000000.0, tolerance);
        // KRATOS_CHECK_RELATIVE_NEAR(rA(5,5), 1.0, tolerance);
    }

}  // namespace Kratos::Testing.