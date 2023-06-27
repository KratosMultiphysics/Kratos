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
#include "trilinos_cpp_test_utilities.h"
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
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    // Linear solvers Trilinos definitions
    using TrilinosLinearSolverType = LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using AmgclMPISolverType = AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

    // Builder and solvers definition
    using TrilinosBuilderAndSolverType = BuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    using TrilinosBlockBuilderAndSolverType = TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    using TrilinosResidualBasedEliminationBuilderAndSolverType = TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;

    // The time scheme
    using TrilinosSchemeType = Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using TrilinosResidualBasedIncrementalUpdateStaticSchemeType = ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void BasicTestBuilderAndSolverDisplacement(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator,
        const bool WithConstraint = false,
        const bool AdditionalNode = false
        )
    {
        // Set MPI coomunicator
        ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, rDataCommunicator);

        // Add variables
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
        if (world_size == 1) {
            // Create nodes
            auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            auto pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            // Create elements
            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
        } else {
            if (rank == 0) {
                // Create nodes
                auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
                pnode1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            } else if (rank == 1) {
                // Create nodes
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                auto pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

                /// Add PARTITION_INDEX
                for (auto& r_node : rModelPart.Nodes()) {
                    r_node.FastGetSolutionStepValue(PARTITION_INDEX) = 1;
                }

                // Create elements
                GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
            }
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

        // TODO: Fix this once MPC work!
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

        // Compute communicaton plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();
    }

    /**
    * @brief It generates a truss structure with an expected solution with an element with null contribution
    */
    static inline void BasicTestBuilderAndSolverDisplacementWithZeroContribution(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Set MPI coomunicator
        ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, rDataCommunicator);

        // Add variables
        rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);
        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(REACTION);
        rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Create properties
        auto p_prop_1 = rModelPart.CreateNewProperties(1, 0);
        p_prop_1->SetValue(YOUNG_MODULUS, 206900000000.0);
        p_prop_1->SetValue(NODAL_AREA, 0.01);

        auto p_prop_2 = rModelPart.CreateNewProperties(2, 0);
        p_prop_2->SetValue(YOUNG_MODULUS, 0.0);
        p_prop_2->SetValue(NODAL_AREA, 0.0);

        // MPI data
        const DataCommunicator& r_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int rank =  r_comm.Rank();
        const int world_size = r_comm.Size();

        // Initially everything in one partition
        if (world_size == 1) {
            NodeType::Pointer pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop_1));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop_2));
        } else {
            if (rank == 0) {
                // Create nodes
                auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
                pnode1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop_1));
            } else if (rank == 1) {
                // Create nodes
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
                auto pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
                pnode3->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop_2));
            }
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

        // Compute communicaton plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();
    }

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void ExtendedTestBuilderAndSolverDisplacement(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator, 
        const bool WithConstraint = false
        )
    {
        // Set MPI coomunicator
        ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, rDataCommunicator);

        // Add variables
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
        if (world_size == 1) {
            NodeType::Pointer pnode1 = rModelPart.CreateNewNode(1, 10.0, -5.0, 0.0);
            NodeType::Pointer pnode2 = rModelPart.CreateNewNode(2, 8.0, -4.0, 0.0);
            NodeType::Pointer pnode3 = rModelPart.CreateNewNode(3, 6.0, -3.0, 0.0);
            NodeType::Pointer pnode4 = rModelPart.CreateNewNode(4, 10.0, 0.0, 0.0);
            NodeType::Pointer pnode5 = rModelPart.CreateNewNode(5, 8.0, 0.0, 0.0);
            NodeType::Pointer pnode6 = rModelPart.CreateNewNode(6, 6.0, 0.0, 0.0);
            NodeType::Pointer pnode7 = rModelPart.CreateNewNode(7, 4.0, -2.0, 0.0);
            NodeType::Pointer pnode8 = rModelPart.CreateNewNode(8, 4.0, 0.0, 0.0);
            NodeType::Pointer pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
            NodeType::Pointer pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);
            NodeType::Pointer pnode11 = rModelPart.CreateNewNode(11, 0.0, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode11, pnode10})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode10, pnode8})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
            GeometryType::Pointer pgeom3 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode8, pnode6})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 3, pgeom3, p_prop));
            GeometryType::Pointer pgeom4 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode6, pnode5})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 4, pgeom4, p_prop));
            GeometryType::Pointer pgeom5 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode5, pnode4})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 5, pgeom5, p_prop));
            GeometryType::Pointer pgeom6 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode4, pnode1})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 6, pgeom6, p_prop));
            GeometryType::Pointer pgeom7 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 7, pgeom7, p_prop));
            GeometryType::Pointer pgeom8 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 8, pgeom8, p_prop));
            GeometryType::Pointer pgeom9 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode3, pnode7})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 9, pgeom9, p_prop));
            GeometryType::Pointer pgeom10 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode7, pnode9})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 10, pgeom10, p_prop));
            GeometryType::Pointer pgeom11 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode9, pnode11})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 11, pgeom11, p_prop));
            GeometryType::Pointer pgeom12 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode10, pnode9})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 12, pgeom12, p_prop));
            GeometryType::Pointer pgeom13 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode9, pnode8})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 13, pgeom13, p_prop));
            GeometryType::Pointer pgeom14 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode8, pnode7})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 14, pgeom14, p_prop));
            GeometryType::Pointer pgeom15 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode7, pnode6})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 15, pgeom15, p_prop));
            GeometryType::Pointer pgeom16 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode6, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 16, pgeom16, p_prop));
            GeometryType::Pointer pgeom17 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode3, pnode5})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 17, pgeom17, p_prop));
            GeometryType::Pointer pgeom18 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode5, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 18, pgeom18, p_prop));
            GeometryType::Pointer pgeom19 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode4})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 19, pgeom19, p_prop));
        } else { // if (world_size == 1) { // TODO: Do more than one partition
            if (rank == 0) {
                NodeType::Pointer pnode1 = rModelPart.CreateNewNode(1, 10.0, -5.0, 0.0);
                NodeType::Pointer pnode2 = rModelPart.CreateNewNode(2, 8.0, -4.0, 0.0);
                NodeType::Pointer pnode3 = rModelPart.CreateNewNode(3, 6.0, -3.0, 0.0);
                NodeType::Pointer pnode4 = rModelPart.CreateNewNode(4, 10.0, 0.0, 0.0);
                NodeType::Pointer pnode5 = rModelPart.CreateNewNode(5, 8.0, 0.0, 0.0);
                NodeType::Pointer pnode6 = rModelPart.CreateNewNode(6, 6.0, 0.0, 0.0);
                NodeType::Pointer pnode7 = rModelPart.CreateNewNode(7, 4.0, -2.0, 0.0);
                NodeType::Pointer pnode8 = rModelPart.CreateNewNode(8, 4.0, 0.0, 0.0);
                NodeType::Pointer pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
                NodeType::Pointer pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);

                /// Add PARTITION_INDEX
                for (auto& r_node : rModelPart.Nodes()) {
                    r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
                }
                pnode9->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
                pnode10->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode10, pnode8})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
                GeometryType::Pointer pgeom3 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode8, pnode6})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 3, pgeom3, p_prop));
                GeometryType::Pointer pgeom4 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode6, pnode5})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 4, pgeom4, p_prop));
                GeometryType::Pointer pgeom5 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode5, pnode4})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 5, pgeom5, p_prop));
                GeometryType::Pointer pgeom6 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode4, pnode1})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 6, pgeom6, p_prop));
                GeometryType::Pointer pgeom7 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode1, pnode2})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 7, pgeom7, p_prop));
                GeometryType::Pointer pgeom8 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode3})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 8, pgeom8, p_prop));
                GeometryType::Pointer pgeom9 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode3, pnode7})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 9, pgeom9, p_prop));
                GeometryType::Pointer pgeom10 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode7, pnode9})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 10, pgeom10, p_prop));
                GeometryType::Pointer pgeom12 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode10, pnode9})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 12, pgeom12, p_prop));
                GeometryType::Pointer pgeom13 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode9, pnode8})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 13, pgeom13, p_prop));
                GeometryType::Pointer pgeom14 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode8, pnode7})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 14, pgeom14, p_prop));
                GeometryType::Pointer pgeom15 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode7, pnode6})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 15, pgeom15, p_prop));
                GeometryType::Pointer pgeom16 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode6, pnode3})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 16, pgeom16, p_prop));
                GeometryType::Pointer pgeom17 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode3, pnode5})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 17, pgeom17, p_prop));
                GeometryType::Pointer pgeom18 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode5, pnode2})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 18, pgeom18, p_prop));
                GeometryType::Pointer pgeom19 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode2, pnode4})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 19, pgeom19, p_prop));
            } else if (rank == 1) {
                NodeType::Pointer pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
                NodeType::Pointer pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);
                NodeType::Pointer pnode11 = rModelPart.CreateNewNode(11, 0.0, 0.0, 0.0);

                /// Add PARTITION_INDEX
                for (auto& r_node : rModelPart.Nodes()) {
                    r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
                }

                GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode11, pnode10})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
                GeometryType::Pointer pgeom11 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({pnode9, pnode11})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 11, pgeom11, p_prop));
            }
        }

        /// Add dof
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X, REACTION_X);
            r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        }

        /// Initialize elements
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& elem : rModelPart.Elements()) {
            elem.Initialize(r_process_info);
            elem.InitializeSolutionStep(r_process_info);
        }

        // Set initial solution
        for (auto& r_node : rModelPart.Nodes()) {
            (r_node.FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (r_node.FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
            (r_node.FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();
        }

        // Fix dofs
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.Fix(DISPLACEMENT_Z);
        }
        // Fix X in first node
        if (rModelPart.HasNode(1)) {
            auto pnode = rModelPart.pGetNode(1);
            pnode->Fix(DISPLACEMENT_X);
        }
        // Fix X and Y in fourth node
        if (rModelPart.HasNode(4)) {
            auto pnode = rModelPart.pGetNode(4);
            pnode->Fix(DISPLACEMENT_X);
            pnode->Fix(DISPLACEMENT_Y);
        }

        // // TODO: Fix this once MPC work!
        // if (WithConstraint) {
        //     rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode1, DISPLACEMENT_Y, *pnode2, DISPLACEMENT_Y, 1.0, 0.0);
        // }

        // Compute communicaton plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();
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

    // static void DebugLHS(const TrilinosSparseSpaceType::MatrixType& rA)
    // {
    //     std::vector<int> row_indexes;
    //     std::vector<int> column_indexes;
    //     std::vector<double> values;
    //     TrilinosCPPTestUtilities::GenerateSparseMatrixIndexAndValuesVectors(rA, row_indexes, column_indexes, values, true, 0.99);
    // }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
    {
        // The base model part
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

        // The data communicator
        const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        // Basic build
        BasicTestBuilderAndSolverDisplacement(r_model_part, r_comm);

        // Create the solvers and things required
        auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
        auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
        Parameters parameters = Parameters(R"(
        {
            "diagonal_values_for_dirichlet_dofs" : "no_scaling",
            "guess_row_size"                     : 15,
            "silent_warnings"                    : false
        })" );
        auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

        // // To create the solution of reference
        // DebugLHS(rA);

        // The solution check
        KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 28);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 2, 3, 4, 4, 5};
        std::vector<int> column_indexes = {0, 1, 2, 4, 3, 2, 4, 5};
        std::vector<double> values = {2069000000.0, 1.0, 4138000000.0, -2069000000.0, 1.0, -2069000000.0, 2069000000.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Testing scale
        parameters["diagonal_values_for_dirichlet_dofs"].SetString("defined_in_process_info");
        r_model_part.GetProcessInfo().SetValue(BUILD_SCALE_FACTOR, 2.26648e+10);
        auto p_builder_and_solver_scale = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA_scale = BuildSystem(r_model_part, p_scheme, p_builder_and_solver_scale);

        // // To create the solution of reference
        // DebugLHS(rA_scale, r_comm);

        KRATOS_CHECK_EQUAL(rA_scale.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA_scale.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA_scale.NumGlobalNonzeros(), 28);

        // Values to check
        row_indexes = {0, 1, 2, 2, 3, 4, 4, 5};
        column_indexes = {0, 1, 2, 4, 3, 2, 4, 5};
        values = {2069000000.0, 22664800000.0, 4138000000.0, -2069000000.0, 22664800000.0, -2069000000.0, 2069000000.0, 22664800000.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA_scale, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system with zero contribution
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverWithZeroContribution, KratosTrilinosApplicationMPITestSuite)
    {
        // The base model part
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

        // The data communicator
        const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        // Basic build
        BasicTestBuilderAndSolverDisplacementWithZeroContribution(r_model_part, r_comm);

        // Create the solvers and things required
        auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
        auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
        Parameters parameters = Parameters(R"(
        {
            "diagonal_values_for_dirichlet_dofs" : "no_scaling",
            "guess_row_size"                     : 15,
            "silent_warnings"                    : false
        })" );
        auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

        // // To create the solution of reference
        // DebugLHS(rA);

        // The solution check
        KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 28);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5};
        std::vector<double> values = {2069000000.0, 1.0, 2069000000.0, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Testing scale
        parameters["diagonal_values_for_dirichlet_dofs"].SetString("defined_in_process_info");
        r_model_part.GetProcessInfo().SetValue(BUILD_SCALE_FACTOR, 2.26648e+10);
        auto p_builder_and_solver_scale = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA_scale = BuildSystem(r_model_part, p_scheme, p_builder_and_solver_scale);

        // // To create the solution of reference
        // DebugLHS(rA_scale, r_comm);

        KRATOS_CHECK_EQUAL(rA_scale.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA_scale.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA_scale.NumGlobalNonzeros(), 28);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5};
        column_indexes = {0, 1, 2, 3, 4, 5};
        values = {2069000000.0, 22664800000.0, 2069000000.0, 22664800000.0, 22664800000.0, 22664800000.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA_scale, row_indexes, column_indexes, values);
    }

    // NOTE: Fails with more than one partition
    // /**
    // * Checks if the elimination builder and solver performs correctly the assemble of the system
    // */
    // KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementEliminationBuilderAndSolverWithZeroContribution, KratosTrilinosApplicationMPITestSuite)
    // {
    //     // The base model part
    //     Model current_model;
    //     ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

    //     // The data communicator
    //     const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

    //     // Generate Epetra communicator
    //     KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    //     auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    //     Epetra_MpiComm epetra_comm(raw_mpi_comm);

    //     // Basic build
    //     BasicTestBuilderAndSolverDisplacement(r_model_part, r_comm);

    //     // Create the solvers and things required
    //     auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
    //     auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
    //     auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosResidualBasedEliminationBuilderAndSolverType(epetra_comm, 15, p_solver) );

    //     const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

    //     // // To create the solution of reference
    //     // DebugLHS(rA);

    //     // The solution check
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 2);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 2);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 4);

    //     // Values to check
    //     std::vector<int> row_indexes = {0, 0, 1, 1};
    //     std::vector<int> column_indexes = {0, 1, 0, 1};
    //     std::vector<double> values = {4138000000.0, -2069000000.0, -2069000000.0, 2069000000.0};

    //     // Check assembly
    //     TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
    // }

    // NOTE: Fails with more than one partition
    // /**
    // * Checks if the elimination builder and solver performs correctly the assemble of the system with zero contribution
    // */
    // KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementEliminationBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
    // {
    //     // The base model part
    //     Model current_model;
    //     ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

    //     // The data communicator
    //     const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

    //     // Generate Epetra communicator
    //     KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    //     auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    //     Epetra_MpiComm epetra_comm(raw_mpi_comm);

    //     // Basic build
    //     BasicTestBuilderAndSolverDisplacementWithZeroContribution(r_model_part, r_comm);

    //     // Create the solvers and things required
    //     auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
    //     auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
    //     auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosResidualBasedEliminationBuilderAndSolverType(epetra_comm, 15, p_solver) );

    //     const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

    //     // // To create the solution of reference
    //     // DebugLHS(rA);

    //     // The solution check
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 2);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 2);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 4);

    //     // Values to check
    //     std::vector<int> row_indexes = {0};
    //     std::vector<int> column_indexes = {0};
    //     std::vector<double> values = {2069000000.0};

    //     // Check assembly
    //     TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
    // }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the extended system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosExtendedDisplacementBlockBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
    {
        // The base model part
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

        // The data communicator
        const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        // Extended build
        ExtendedTestBuilderAndSolverDisplacement(r_model_part, r_comm);

        // Create the solvers and things required
        auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
        auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
        Parameters parameters = Parameters(R"(
        {
            "diagonal_values_for_dirichlet_dofs" : "no_scaling",
            "guess_row_size"                     : 15,
            "silent_warnings"                    : false
        })" );
        auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

        // // To create the solution of reference
        // DebugLHS(rA);

        // The solution check
        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 22);
        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 196);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21};
        std::vector<int> column_indexes = {0, 1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 9, 2, 3, 4, 5, 8, 9, 12, 13, 2, 3, 4, 5, 8, 9, 11, 12, 13, 6, 7, 4, 5, 8, 9, 10, 3, 4, 5, 8, 9, 8, 10, 11, 12, 13, 14, 5, 10, 11, 12, 13, 4, 5, 10, 11, 12, 13, 16, 17, 4, 5, 10, 11, 12, 13, 15, 16, 17, 10, 14, 15, 16, 17, 18, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17, 20, 21, 12, 13, 14, 15, 16, 17, 19, 20, 21, 14, 18, 20, 17, 19, 16, 17, 18, 20, 21, 16, 17, 20, 21};
        std::vector<double> values = {740227943.2715302705764771, 598856985.8178827762603760, 370113971.6357653141021729, -185056985.8178827166557312, 370113971.6357653141021729, 1572984379.4520018100738525, -555170957.4536480903625488, -740227943.2715302705764771, 370113971.6357653141021729, -185056985.8178827166557312, -555170957.4536479711532593, 1257477943.2715306282043457, 370113971.6357653141021729, -185056985.8178827166557312, -517250000.0, -740227943.2715302705764771, 370113971.6357653141021729, 1657021225.9261374473571777, -475379934.1969153881072998, -176565339.3830768167972565, -264848009.0746152102947235, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -475379934.1969153881072998, 1457052651.9143548011779785, -264848009.0746152102947235, -397272013.6119228005409241, -689666666.6666666269302368, 370113971.6357653141021729, -185056985.8178827166557312, 1127028492.9089412689208984, 783913971.6357650756835938, -176565339.3830768167972565, -264848009.0746152102947235, 2245565339.3830766677856445, 264848009.0746151804924011, -1034500000.0, -517250000.0, -264848009.0746152102947235, -397272013.6119228005409241, 264848009.0746151804924011, 914522013.6119227409362793, -1034500000.0, 2434750982.5687417984008789, 365750982.5687416195869446, -365750982.5687417387962341, -365750982.5687416791915894, -1034500000.0, -689666666.6666666269302368, 365750982.5687416195869446, 1055417649.2354083061218262, -365750982.5687416791915894, -365750982.5687416195869446, -740227943.2715302705764771, 370113971.6357653141021729, -365750982.5687417387962341, -365750982.5687416791915894, 1846206869.1118023395538330, -374476960.7027890682220459, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -365750982.5687416791915894, -365750982.5687416195869446, -374476960.7027890682220459, 1770364954.2045071125030518, -1034500000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2809227943.2715301513671875, 370113971.6357650756835938, -740227943.2715302705764771, -370113971.6357651352882385, -1034500000.0, -1034500000.0, 370113971.6357650756835938, 1219556985.8178825378417969, -370113971.6357651352882385, -185056985.8178825676441193, -740227943.2715302705764771, 370113971.6357653141021729, -740227943.2715302705764771, -370113971.6357651352882385, 2220683829.8145909309387207, -370113971.6357656121253967, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357651352882385, -185056985.8178825676441193, -370113971.6357656121253967, 2624170957.4536480903625488, -2069000000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2069000000.0, -1034500000.0, -2069000000.0, 2069000000.0, -740227943.2715302705764771, 370113971.6357653141021729, -1034500000.0, 1774727943.2715301513671875, -370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357653141021729, 185056985.8178827166557312};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
    }

    // NOTE: Fails with more than one partition
    // /**
    // * Checks if the elimination builder and solver performs correctly the assemble of the extended system
    // */
    // KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosExtendedDisplacementEliminationBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
    // {
    //     // The base model part
    //     Model current_model;
    //     ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

    //     // The data communicator
    //     const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

    //     // Generate Epetra communicator
    //     KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    //     auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
    //     Epetra_MpiComm epetra_comm(raw_mpi_comm);

    //     // Extended build
    //     ExtendedTestBuilderAndSolverDisplacement(r_model_part, r_comm);

    //     // Create the solvers and things required
    //     auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
    //     auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
    //     auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosResidualBasedEliminationBuilderAndSolverType(epetra_comm, 15, p_solver) );

    //     const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

    //     // // To create the solution of reference
    //     // DebugLHS(rA);

    //     // The solution check
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 19);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 19);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 161);

    //     // Values to check
    //     std::vector<int> row_indexes = {0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18};
    //     std::vector<int> column_indexes = {0, 1, 2, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 6, 1, 2, 3, 4, 5, 6, 9, 10, 1, 2, 3, 4, 5, 6, 8, 9, 10, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 5, 7, 8, 9, 10, 11, 4, 7, 8, 9, 10, 3, 4, 7, 8, 9, 10, 13, 14, 3, 4, 7, 8, 9, 10, 12, 13, 14, 7, 11, 12, 13, 14, 15, 10, 11, 12, 13, 14, 9, 10, 11, 12, 13, 14, 17, 18, 9, 10, 11, 12, 13, 14, 16, 17, 18, 11, 15, 17, 14, 16, 13, 14, 15, 17, 18, 13, 14, 17, 18};
    //     std::vector<double> values = {598856985.8178827762603760, 370113971.6357653141021729, -185056985.8178827166557312, 370113971.6357653141021729, 1572984379.4520018100738525, -555170957.4536480903625488, -740227943.2715302705764771, 370113971.6357653141021729, -185056985.8178827166557312, -555170957.4536479711532593, 1257477943.2715306282043457, 370113971.6357653141021729, -185056985.8178827166557312, -517250000.0, -740227943.2715302705764771, 370113971.6357653141021729, 1657021225.9261374473571777, -475379934.1969153881072998, -176565339.3830768167972565, -264848009.0746152102947235, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -475379934.1969153881072998, 1457052651.9143548011779785, -264848009.0746152102947235, -397272013.6119228005409241, -689666666.6666666269302368, 370113971.6357653141021729, -185056985.8178827166557312, -176565339.3830768167972565, -264848009.0746152102947235, 2245565339.3830766677856445, 264848009.0746151804924011, -1034500000.0, -517250000.0, -264848009.0746152102947235, -397272013.6119228005409241, 264848009.0746151804924011, 914522013.6119227409362793, -1034500000.0, 2434750982.5687417984008789, 365750982.5687416195869446, -365750982.5687417387962341, -365750982.5687416791915894, -1034500000.0, -689666666.6666666269302368, 365750982.5687416195869446, 1055417649.2354083061218262, -365750982.5687416791915894, -365750982.5687416195869446, -740227943.2715302705764771, 370113971.6357653141021729, -365750982.5687417387962341, -365750982.5687416791915894, 1846206869.1118023395538330, -374476960.7027890682220459, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -365750982.5687416791915894, -365750982.5687416195869446, -374476960.7027890682220459, 1770364954.2045071125030518, -1034500000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2809227943.2715301513671875, 370113971.6357650756835938, -740227943.2715302705764771, -370113971.6357651352882385, -1034500000.0, -1034500000.0, 370113971.6357650756835938, 1219556985.8178825378417969, -370113971.6357651352882385, -185056985.8178825676441193, -740227943.2715302705764771, 370113971.6357653141021729, -740227943.2715302705764771, -370113971.6357651352882385, 2220683829.8145909309387207, -370113971.6357656121253967, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357651352882385, -185056985.8178825676441193, -370113971.6357656121253967, 2624170957.4536480903625488, -2069000000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2069000000.0, -1034500000.0, -2069000000.0, 2069000000.0, -740227943.2715302705764771, 370113971.6357653141021729, -1034500000.0, 1774727943.2715301513671875, -370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357653141021729, 185056985.8178827166557312};

    //     // Check assembly
    //     TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
    // }

}  // namespace Kratos::Testing.