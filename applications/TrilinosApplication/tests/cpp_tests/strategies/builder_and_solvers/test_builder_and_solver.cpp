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
#include "../../trilinos_cpp_test_utilities.h"
#include "containers/model.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "mpi/tests/test_utilities/mpi_cpp_test_utilities.h"

/* Element include */
#include "geometries/line_2d_2.h"
#include "tests/test_utilities/test_bar_element.h"

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
    using GeometryType = Geometry<Node>;
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    /// Linear solvers Trilinos definitions
    using TrilinosLinearSolverType = LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using AmgclMPISolverType = AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

    /// Builder and solvers definition
    using TrilinosBuilderAndSolverType = BuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    using TrilinosBlockBuilderAndSolverType = TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    using TrilinosResidualBasedEliminationBuilderAndSolverType = TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;

    /// The time scheme
    using TrilinosSchemeType = Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    using TrilinosResidualBasedIncrementalUpdateStaticSchemeType = ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void DynamicAllocationBasicTestBuilderAndSolverDisplacement(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator,
        const bool WithConstraint = false
        )
    {
        // Set MPI communicator
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
        const int rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

        // Create two nodes per partition
        Node::Pointer pnode0 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        Node::Pointer pnode1 = rank == 0 ? pnode0 : rModelPart.CreateNewNode(rank + 1, static_cast<double>(rank), 0.0, 0.0);
        Node::Pointer pnode2 = rModelPart.CreateNewNode(rank + 2, static_cast<double>(rank + 1), 0.0, 0.0);

        /// Add PARTITION_INDEX
        pnode0->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
        pnode1->FastGetSolutionStepValue(PARTITION_INDEX) = rank;
        if (rank + 1 < world_size) {
            pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = rank + 1;
        } else {
            pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = rank;
        }

        // Create elements
        GeometryType::Pointer pgeom = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( rank + 1, pgeom, p_prop));


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

        // Adding constraints
        if (WithConstraint) {
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", (2 * rank) + 1, *pnode0, DISPLACEMENT_X, *pnode2, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", (2 * rank) + 2, *pnode0, DISPLACEMENT_Y, *pnode2, DISPLACEMENT_Y, 1.0, 0.0);
        }

        // Compute communication plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();
    }

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void BasicTestBuilderAndSolverDisplacement(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator,
        const bool WithConstraint = false,
        const bool AdditionalNode = false,
        const bool InvertRoleAdditionalNode = false
        )
    {
        // Set MPI communicator
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
        const int rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

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
            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
        } else {
            if (rank == 0) {
                // Create nodes
                auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
                pnode1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
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
                GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
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

        // Adding constraints
        if (WithConstraint) {
            if (rModelPart.HasNode(3)) {
                auto pnode3 = rModelPart.pGetNode(3);
                const auto node_3_partition_index = pnode3->FastGetSolutionStepValue(PARTITION_INDEX);
                if (rModelPart.HasNode(2)) {
                    auto pnode2 = rModelPart.pGetNode(2);
                    rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode2, DISPLACEMENT_X, *pnode3, DISPLACEMENT_X, 1.0, 0.0);
                }

                if (AdditionalNode) {
                    auto pnode4 = rModelPart.CreateNewNode(4, 3.0, 0.0, 0.0);
                    pnode4->FastGetSolutionStepValue(PARTITION_INDEX) = node_3_partition_index;
                    pnode4->AddDof(DISPLACEMENT_X, REACTION_X);
                    pnode4->AddDof(DISPLACEMENT_Y, REACTION_Y);
                    pnode4->AddDof(DISPLACEMENT_Z, REACTION_Z);
                    auto pnode5 = rModelPart.CreateNewNode(5, 4.0, 0.0, 0.0);
                    pnode5->FastGetSolutionStepValue(PARTITION_INDEX) = node_3_partition_index;
                    pnode5->AddDof(DISPLACEMENT_X, REACTION_X);
                    pnode5->AddDof(DISPLACEMENT_Y, REACTION_Y);
                    pnode5->AddDof(DISPLACEMENT_Z, REACTION_Z);
                    // Existing node is the master
                    if (!InvertRoleAdditionalNode) {
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, *pnode3, DISPLACEMENT_X, *pnode4, DISPLACEMENT_X, 1.0, 0.0);
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, *pnode3, DISPLACEMENT_X, *pnode5, DISPLACEMENT_X, 1.0, 0.0);
                    } else { // Free nodes are the master
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, *pnode4, DISPLACEMENT_X, *pnode3, DISPLACEMENT_X, 1.0, 0.0);
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, *pnode5, DISPLACEMENT_X, *pnode3, DISPLACEMENT_X, 1.0, 0.0);
                    }
                }
            }
        }

        // Compute communication plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();
    }

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void BasicTestBuilderAndSolverDisplacementAllDoFsMaster(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator,
        const bool MasterElementBelongsToStructure = false
        )
    {
        // Set MPI communicator
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
        const int rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

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
            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
        } else {
            if (rank == 0) {
                // Create nodes
                auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
                pnode1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
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
                GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
            }
        }

        // Adding fourth node
        Node::Pointer pnode4 = rModelPart.CreateNewNode(4, 3.0, 0.0, 0.0);
        const int rank_point_4 = (world_size == 1) ? 0 : 1;
        pnode4->FastGetSolutionStepValue(PARTITION_INDEX) = rank_point_4;

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
        if (MasterElementBelongsToStructure) {
            if(rModelPart.HasNode(2)) {
                auto pnode2 = rModelPart.pGetNode(2);
                pnode2->Fix(DISPLACEMENT_X);
                pnode2->Fix(DISPLACEMENT_Y);
                pnode2->Fix(DISPLACEMENT_Z);
                for (auto& r_node : rModelPart.Nodes()) {
                    if (r_node.Id() != 2) {
                        const auto i = r_node.Id() - 1;
                        if (rank == r_node.FastGetSolutionStepValue(PARTITION_INDEX)) {
                            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3 * i + 1, *pnode2, DISPLACEMENT_X, r_node, DISPLACEMENT_X, 1.0, 0.0);
                            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3 * i + 2, *pnode2, DISPLACEMENT_Y, r_node, DISPLACEMENT_Y, 1.0, 0.0);
                            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3 * i + 3, *pnode2, DISPLACEMENT_Z, r_node, DISPLACEMENT_Z, 1.0, 0.0);
                        }
                    }
                }
            }
        } else {
            pnode4->Fix(DISPLACEMENT_X);
            pnode4->Fix(DISPLACEMENT_Y);
            pnode4->Fix(DISPLACEMENT_Z);
            for (auto& r_node : rModelPart.Nodes()) {
                if (r_node.Id() != 4) {
                    const auto i = r_node.Id() - 1;
                    if (rank == r_node.FastGetSolutionStepValue(PARTITION_INDEX)) {
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3 * i + 1, *pnode4, DISPLACEMENT_X, r_node, DISPLACEMENT_X, 1.0, 0.0);
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3 * i + 2, *pnode4, DISPLACEMENT_Y, r_node, DISPLACEMENT_Y, 1.0, 0.0);
                        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3 * i + 3, *pnode4, DISPLACEMENT_Z, r_node, DISPLACEMENT_Z, 1.0, 0.0);
                    }
                }
            }
        }

        // Compute communication plan and fill communicator meshes correctly
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
        // Set MPI communicator
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
        const int rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

        // Initially everything in one partition
        if (world_size == 1) {
            Node::Pointer pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            Node::Pointer pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            Node::Pointer pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop_1));
            GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop_2));
        } else {
            if (rank == 0) {
                // Create nodes
                auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
                pnode1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
                rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop_1));
            } else if (rank == 1) {
                // Create nodes
                auto pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
                auto pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
                pnode3->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

                // Create elements
                GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
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

        // Compute communication plan and fill communicator meshes correctly
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
        // Add variables
        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(REACTION);
        rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Fill model part
        MPICppTestUtilities::GenerateDistributedBarStructure(rModelPart, rDataCommunicator);

        // Create properties
        auto p_prop = rModelPart.pGetProperties(1, 0);
        p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
        p_prop->SetValue(NODAL_AREA, 0.01);

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

        // Adding constraint
        if (WithConstraint) {
            if (rModelPart.HasNode(1) && rModelPart.HasNode(2)) {
                auto pnode1 = rModelPart.pGetNode(1);
                auto pnode2 = rModelPart.pGetNode(2);
                rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode1, DISPLACEMENT_Y, *pnode2, DISPLACEMENT_Y, 1.0, 0.0);
            }
        }

        // Compute communication plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();
    }

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline bool StructuralTestBuilderAndSolverDisplacement(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator,
        const bool InvertLocationConstraint = false
        )
    {
        // Check structural app is compiled
        if (!KratosComponents<Element>::Has("SmallDisplacementElement3D8N")) {
            return false;
        }

        // Add variables
        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(REACTION);
        rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);

        // MPI data
        const int rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

        // Fill model part
        const double length_cube = 1.0;
        auto& r_sub_model_part_1 = rModelPart.CreateSubModelPart("1");
        auto& r_sub_model_part_2 = rModelPart.CreateSubModelPart("2");
        if (rank == 0) {
            // Create nodes
            r_sub_model_part_1.CreateNewNode(1, 0.0, length_cube, length_cube);
            r_sub_model_part_1.CreateNewNode(2, 0.0, length_cube, 0.0);
            r_sub_model_part_1.CreateNewNode(3, 0.0, 0.0, length_cube);
            r_sub_model_part_1.CreateNewNode(4, length_cube, length_cube, length_cube);
            r_sub_model_part_1.CreateNewNode(5, 0.0, 0.0, 0.0);
            r_sub_model_part_1.CreateNewNode(6, length_cube, length_cube, 0.0);
            r_sub_model_part_1.CreateNewNode(7, length_cube, 0.0, length_cube);
            r_sub_model_part_1.CreateNewNode(8, length_cube, 0.0, 0.0);
            r_sub_model_part_1.CreateNewNode(9, 2 * length_cube, length_cube, length_cube);
            r_sub_model_part_1.CreateNewNode(10, 2 * length_cube, length_cube, 0.0);
            r_sub_model_part_1.CreateNewNode(11, 2 * length_cube, 0.0, length_cube);
            r_sub_model_part_1.CreateNewNode(12, 2 * length_cube, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            // Other nodes
            if (!InvertLocationConstraint || world_size == 1) {
                const int other_nodes_rank = world_size == 1 ? 0 : 1;
                r_sub_model_part_2.CreateNewNode(13, 0.0, length_cube, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_2.CreateNewNode(14, 0.0, 0.0, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_2.CreateNewNode(15, length_cube, length_cube, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_2.CreateNewNode(16, length_cube, 0.0, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_2.CreateNewNode(17, 2 * length_cube, length_cube, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_2.CreateNewNode(18, 2 * length_cube, 0.0, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
            }
        } else if (rank == 1) {
            // Create nodes
            r_sub_model_part_2.CreateNewNode(13, 0.0, length_cube, length_cube);
            r_sub_model_part_2.CreateNewNode(14, 0.0, 0.0, length_cube);
            r_sub_model_part_2.CreateNewNode(15, length_cube, length_cube, length_cube);
            r_sub_model_part_2.CreateNewNode(16, length_cube, 0.0, length_cube);
            r_sub_model_part_2.CreateNewNode(17, 2 * length_cube, length_cube, length_cube);
            r_sub_model_part_2.CreateNewNode(18, 2 * length_cube, 0.0, length_cube);

            /// Add PARTITION_INDEX
            for (auto& r_node : r_sub_model_part_2.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            // Other nodes
            if (InvertLocationConstraint) {
                const int other_nodes_rank = 0;
                r_sub_model_part_1.CreateNewNode(1, 0.0, length_cube, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_1.CreateNewNode(3, 0.0, 0.0, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_1.CreateNewNode(4, length_cube, length_cube, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_1.CreateNewNode(7, length_cube, 0.0, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_1.CreateNewNode(9, 2 * length_cube, length_cube, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
                r_sub_model_part_1.CreateNewNode(11, 2 * length_cube, 0.0, length_cube)->FastGetSolutionStepValue(PARTITION_INDEX) = other_nodes_rank;
            }
        }

        /// Add dof
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X, REACTION_X);
            r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        }

        // Add elements
        auto p_prop0 = rModelPart.CreateNewProperties(0);
        auto p_prop1 = rModelPart.CreateNewProperties(1);
        if (rank == 0) {
            // Create properties
            p_prop0->SetValue(YOUNG_MODULUS,1.0e6);
            p_prop0->SetValue(POISSON_RATIO,0.0);
            p_prop0->SetValue(DENSITY,1.0);
            const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            p_prop0->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            // Create elements
            rModelPart.CreateNewElement("SmallDisplacementElement3D8N", 1, {5,8,6,2,3,7,4,1}, p_prop0);
            rModelPart.CreateNewElement("SmallDisplacementElement3D8N", 2, {8,12,10,6,7,11,9,4}, p_prop0);
        }

        if (rank == 1 || world_size == 1) {
            // Create properties
            p_prop1->SetValue(YOUNG_MODULUS,1.0e6);
            p_prop1->SetValue(POISSON_RATIO,0.0);
            p_prop1->SetValue(DENSITY,1.0);
            p_prop1->SetValue(THICKNESS,0.01);
            const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
            p_prop1->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            // Create elements
            rModelPart.CreateNewElement("MembraneElement3D4N", 3, {14,16,15,13}, p_prop1);
            rModelPart.CreateNewElement("MembraneElement3D4N", 4, {16,18,17,15}, p_prop1);
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
        for (auto& r_node : r_sub_model_part_1.Nodes()) {
            if (r_node.X() < 0.01) {
                r_node.Fix(DISPLACEMENT_X);
                r_node.Fix(DISPLACEMENT_Y);
                r_node.Fix(DISPLACEMENT_Z);
            } else if (r_node.X() > 2 * length_cube - 0.01) {
                r_node.Fix(DISPLACEMENT_X);
                r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = 0.01;
                r_node.Fix(DISPLACEMENT_Y);
                r_node.Fix(DISPLACEMENT_Z);
            }
        }

        // Adding constraint
        bool add_constraint = false;
        if (InvertLocationConstraint) {
            add_constraint = (rank == 1) || (world_size == 1);
        } else {
            add_constraint = (rank == 0) || (world_size == 1);
        }
        if (add_constraint) {
            auto pnode1 = rModelPart.pGetNode(1);
            auto pnode13 = rModelPart.pGetNode(13);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode1, DISPLACEMENT_X, *pnode13, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, *pnode1, DISPLACEMENT_Y, *pnode13, DISPLACEMENT_Y, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 3, *pnode1, DISPLACEMENT_Z, *pnode13, DISPLACEMENT_Z, 1.0, 0.0);
            auto pnode3 = rModelPart.pGetNode(3);
            auto pnode14 = rModelPart.pGetNode(14);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 4, *pnode3, DISPLACEMENT_X, *pnode14, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 5, *pnode3, DISPLACEMENT_Y, *pnode14, DISPLACEMENT_Y, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 6, *pnode3, DISPLACEMENT_Z, *pnode14, DISPLACEMENT_Z, 1.0, 0.0);
            auto pnode4 = rModelPart.pGetNode(4);
            auto pnode15 = rModelPart.pGetNode(15);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 7, *pnode4, DISPLACEMENT_X, *pnode15, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 8, *pnode4, DISPLACEMENT_Y, *pnode15, DISPLACEMENT_Y, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 9, *pnode4, DISPLACEMENT_Z, *pnode15, DISPLACEMENT_Z, 1.0, 0.0);
            auto pnode7 = rModelPart.pGetNode(7);
            auto pnode16 = rModelPart.pGetNode(16);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 10, *pnode7, DISPLACEMENT_X, *pnode16, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 11, *pnode7, DISPLACEMENT_Y, *pnode16, DISPLACEMENT_Y, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 12, *pnode7, DISPLACEMENT_Z, *pnode16, DISPLACEMENT_Z, 1.0, 0.0);
            auto pnode9 = rModelPart.pGetNode(9);
            auto pnode17 = rModelPart.pGetNode(17);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 13, *pnode9, DISPLACEMENT_X, *pnode17, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 14, *pnode9, DISPLACEMENT_Y, *pnode17, DISPLACEMENT_Y, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 15, *pnode9, DISPLACEMENT_Z, *pnode17, DISPLACEMENT_Z, 1.0, 0.0);
            auto pnode11 = rModelPart.pGetNode(11);
            auto pnode18 = rModelPart.pGetNode(18);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 16, *pnode11, DISPLACEMENT_X, *pnode18, DISPLACEMENT_X, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 17, *pnode11, DISPLACEMENT_Y, *pnode18, DISPLACEMENT_Y, 1.0, 0.0);
            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 18, *pnode11, DISPLACEMENT_Z, *pnode18, DISPLACEMENT_Z, 1.0, 0.0);
        }

        // Compute communication plan and fill communicator meshes correctly
        ParallelFillCommunicator(rModelPart, rDataCommunicator).Execute();

        return true;
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

        TrilinosSparseSpaceType::VectorPointerType pDx; /// The increment in the solution
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
        if(rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            pBuilderAndSolver->ApplyConstraints(pScheme, rModelPart, rA, rb);
        }
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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 28);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 2, 3, 4, 4, 5};
        std::vector<int> column_indexes = {0, 1, 2, 4, 3, 2, 4, 5};
        std::vector<double> values = {E_100, 1.0, 2 * E_100, -E_100, 1.0, -E_100, E_100, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Testing scale
        parameters["diagonal_values_for_dirichlet_dofs"].SetString("defined_in_process_info");
        r_model_part.GetProcessInfo().SetValue(BUILD_SCALE_FACTOR, 2.26648e+10);
        auto p_builder_and_solver_scale = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA_scale = BuildSystem(r_model_part, p_scheme, p_builder_and_solver_scale);

        // // To create the solution of reference
        // DebugLHS(rA_scale);

        KRATOS_EXPECT_EQ(rA_scale.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(rA_scale.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(rA_scale.NumGlobalNonzeros(), 28);

        // Values to check
        row_indexes = {0, 1, 2, 2, 3, 4, 4, 5};
        column_indexes = {0, 1, 2, 4, 3, 2, 4, 5};
        values = {E_100, 22664800000.0, 2 * E_100, -E_100, 22664800000.0, -E_100, E_100, 22664800000.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA_scale, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosDynamicAllocationBasicDisplacementBlockBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
    {
        // The base model part
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

        // The data communicator
        const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

        // The world size
        const int world_size = r_comm.Size();

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        // Basic build
        DynamicAllocationBasicTestBuilderAndSolverDisplacement(r_model_part, r_comm);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 2 * (world_size + 1));
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 2 * (world_size + 1));
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 16 + 12 * (world_size - 1));

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 3};
        std::vector<int> column_indexes = {0, 1, 2, 3};
        std::vector<double> values = {E_100, 1.0, E_100, 1.0};
        if (world_size > 1) {
            values[2] *= 2;
        }
        for (int i = 1; i < world_size; ++i) {
            // Diagonal values
            row_indexes.push_back(2 * (world_size + 1) - 2);
            row_indexes.push_back(2 * (world_size + 1) - 1);
            column_indexes.push_back(2 * (world_size + 1) - 2);
            column_indexes.push_back(2 * (world_size + 1) - 1);
            values.push_back(E_100);
            values.push_back(1.0);

            // Non diagonal values
            const unsigned int index_i = 2 * (world_size + 1) - 2;
            const unsigned int index_j = 2 * (world_size + 1) - 4;
            row_indexes.push_back(index_i);
            row_indexes.push_back(index_j);
            column_indexes.push_back(index_j);
            column_indexes.push_back(index_i);
            values.push_back(-E_100);
            values.push_back(-E_100);
        }

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 28);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5};
        std::vector<double> values = {E_100, 1.0, E_100, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Testing scale
        parameters["diagonal_values_for_dirichlet_dofs"].SetString("defined_in_process_info");
        r_model_part.GetProcessInfo().SetValue(BUILD_SCALE_FACTOR, 2.26648e+10);
        auto p_builder_and_solver_scale = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, p_solver, parameters) );

        const auto& rA_scale = BuildSystem(r_model_part, p_scheme, p_builder_and_solver_scale);

        // // To create the solution of reference
        // DebugLHS(rA_scale);

        KRATOS_EXPECT_EQ(rA_scale.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(rA_scale.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(rA_scale.NumGlobalNonzeros(), 28);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5};
        column_indexes = {0, 1, 2, 3, 4, 5};
        values = {E_100, 22664800000.0, E_100, 22664800000.0, 22664800000.0, 22664800000.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA_scale, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver with constraints performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverWithConstraints, KratosTrilinosApplicationMPITestSuite)
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
        BasicTestBuilderAndSolverDisplacement(r_model_part, r_comm, true);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 28);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5};
        std::vector<double> values = {E_100, 1.0, E_100, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 7);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 4, 5};
        column_indexes = {0, 1, 2, 3, 2, 4, 5};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver with constraints performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverWithInactiveConstraints, KratosTrilinosApplicationMPITestSuite)
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
        BasicTestBuilderAndSolverDisplacement(r_model_part, r_comm, true);
        if (r_model_part.HasMasterSlaveConstraint(1)) {
            r_model_part.pGetMasterSlaveConstraint(1)->Set(ACTIVE, false);
        }

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 28);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 2, 3, 4, 4, 5};
        std::vector<int> column_indexes = {0, 1, 2, 4, 3, 2, 4, 5};
        std::vector<double> values = {E_100, 1.0, 2 * E_100, -E_100, 1.0, -E_100, E_100, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 6);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 6);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 7);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 4, 5};
        column_indexes = {0, 1, 2, 3, 2, 4, 5};
        values = {1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosDynamicAllocationBasicDisplacementBlockBuilderAndSolverWithConstraints, KratosTrilinosApplicationMPITestSuite)
    {
        // The base model part
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main", 3);

        // The data communicator
        const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

        // The world size
        const int world_size = r_comm.Size();

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        // Basic build
        DynamicAllocationBasicTestBuilderAndSolverDisplacement(r_model_part, r_comm, true);

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

        // To create the solution of reference
        // DebugLHS(rA);

        // The solution check
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 2 * (world_size + 1));
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 2 * (world_size + 1));
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 16 + 14 * (world_size - 1));

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 3};
        std::vector<int> column_indexes = {0, 1, 2, 3};
        std::vector<double> values = {1.0, 1.0, 1.0, 1.0};
        for (int i = 1; i < world_size; ++i) {
            const unsigned int index_0 = 4 + 2 * (i - 1);
            const unsigned int index_1 = 5 + 2 * (i - 1);
            row_indexes.push_back(index_0);
            row_indexes.push_back(index_1);
            column_indexes.push_back(index_0);
            column_indexes.push_back(index_1);
            values.push_back(1.0);
            values.push_back(1.0);
        }

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // To create the solution of reference
        // DebugLHS(r_T);

        row_indexes = std::vector<int>({0, 1, 2, 2, 3 ,3});
        column_indexes = std::vector<int>({0, 1, 0, 2 , 1, 3});
        values = std::vector<double>({1.0, 1.0, 1.0, 0.0, 1.0, 0.0});
        for (int i = 1; i < world_size; ++i) {
            const unsigned int index_0 = 4 + 2 * (i - 1);
            const unsigned int index_1 = 5 + 2 * (i - 1);
            row_indexes.push_back(index_0);
            row_indexes.push_back(index_0);
            row_indexes.push_back(index_1);
            row_indexes.push_back(index_1);
            column_indexes.push_back(index_0);
            column_indexes.push_back(0);
            column_indexes.push_back(index_1);
            column_indexes.push_back(1);
            values.push_back(0.0);
            values.push_back(1.0);
            values.push_back(0.0);
            values.push_back(1.0);
        }

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver with constraints performs correctly the assemble of the system with auxiliary node
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverWithConstraintsAuxiliarNode, KratosTrilinosApplicationMPITestSuite)
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
        BasicTestBuilderAndSolverDisplacement(r_model_part, r_comm, true, true);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 8);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 8);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 32);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<double> values = {E_100, 1.0, E_100, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 8);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 8);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 11);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5, 6, 7};
        column_indexes = {0, 1, 2, 3, 2, 5, 4, 4};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver with constraints performs correctly the assemble of the system with auxiliary node (inverted)
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverWithConstraintsAuxiliarNodeInverted, KratosTrilinosApplicationMPITestSuite)
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
        BasicTestBuilderAndSolverDisplacement(r_model_part, r_comm, true, true, true);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 8);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 8);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 38);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5, 6, 6, 7, 7};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 6, 7};
        std::vector<double> values = {E_100, 1.0, E_100, 1.0, 1.0, 1.0, E_100, E_100, E_100, E_100};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 8);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 8);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 11);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 4, 4, 5, 6, 7};
        column_indexes = {0, 1, 2, 3, 2, 6, 7, 5, 6, 7};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverAllDoFsMaster, KratosTrilinosApplicationMPITestSuite)
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
        BasicTestBuilderAndSolverDisplacementAllDoFsMaster(r_model_part, r_comm);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 12);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 12);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 43);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        std::vector<double> values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 12);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 12);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 21);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        column_indexes = {9, 10, 11, 9, 10, 11, 9, 10, 11, 9, 10, 11};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the system
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosBasicDisplacementBlockBuilderAndSolverAllDoFsMasterFromStructureSide, KratosTrilinosApplicationMPITestSuite)
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
        BasicTestBuilderAndSolverDisplacementAllDoFsMaster(r_model_part, r_comm, true);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 12);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 12);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 39);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        std::vector<double> values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 12);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 12);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 21);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        column_indexes = {3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
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
    //     KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 2);
    //     KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 2);
    //     KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 4);

    //     // Values to check
    //     auto p_prop = r_model_part.pGetProperties(1);
    //     const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
    //     std::vector<int> row_indexes = {0, 0, 1, 1};
    //     std::vector<int> column_indexes = {0, 1, 0, 1};
    //     std::vector<double> values = {2 * E_100, -E_100, -E_100, E_100};

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
    //     KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 2);
    //     KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 2);
    //     KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 4);

    //     // Values to check
    //     std::vector<int> row_indexes = {0};
    //     std::vector<int> column_indexes = {0};
    //     std::vector<double> values = {E_100};

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
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 22);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 196);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21};
        std::vector<int> column_indexes = {0, 1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 9, 2, 3, 4, 5, 8, 9, 12, 13, 2, 3, 4, 5, 8, 9, 11, 12, 13, 6, 7, 4, 5, 8, 9, 10, 3, 4, 5, 8, 9, 8, 10, 11, 12, 13, 14, 5, 10, 11, 12, 13, 4, 5, 10, 11, 12, 13, 16, 17, 4, 5, 10, 11, 12, 13, 15, 16, 17, 10, 14, 15, 16, 17, 18, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17, 20, 21, 12, 13, 14, 15, 16, 17, 19, 20, 21, 14, 18, 20, 17, 19, 16, 17, 18, 20, 21, 16, 17, 20, 21};
        std::vector<double> values = {740227943.2715302705764771, 598856985.8178827762603760, 370113971.6357653141021729, -185056985.8178827166557312, 370113971.6357653141021729, 1572984379.4520018100738525, -555170957.4536480903625488, -740227943.2715302705764771, 370113971.6357653141021729, -185056985.8178827166557312, -555170957.4536479711532593, 1257477943.2715306282043457, 370113971.6357653141021729, -185056985.8178827166557312, -517250000.0, -740227943.2715302705764771, 370113971.6357653141021729, 1657021225.9261374473571777, -475379934.1969153881072998, -176565339.3830768167972565, -264848009.0746152102947235, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -475379934.1969153881072998, 1457052651.9143548011779785, -264848009.0746152102947235, -397272013.6119228005409241, -689666666.6666666269302368, 370113971.6357653141021729, -185056985.8178827166557312, 1127028492.9089412689208984, 783913971.6357650756835938, -176565339.3830768167972565, -264848009.0746152102947235, 2245565339.3830766677856445, 264848009.0746151804924011, -1034500000.0, -517250000.0, -264848009.0746152102947235, -397272013.6119228005409241, 264848009.0746151804924011, 914522013.6119227409362793, -1034500000.0, 2434750982.5687417984008789, 365750982.5687416195869446, -365750982.5687417387962341, -365750982.5687416791915894, -1034500000.0, -689666666.6666666269302368, 365750982.5687416195869446, 1055417649.2354083061218262, -365750982.5687416791915894, -365750982.5687416195869446, -740227943.2715302705764771, 370113971.6357653141021729, -365750982.5687417387962341, -365750982.5687416791915894, 1846206869.1118023395538330, -374476960.7027890682220459, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -365750982.5687416791915894, -365750982.5687416195869446, -374476960.7027890682220459, 1770364954.2045071125030518, -1034500000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2809227943.2715301513671875, 370113971.6357650756835938, -740227943.2715302705764771, -370113971.6357651352882385, -1034500000.0, -1034500000.0, 370113971.6357650756835938, 1219556985.8178825378417969, -370113971.6357651352882385, -185056985.8178825676441193, -740227943.2715302705764771, 370113971.6357653141021729, -740227943.2715302705764771, -370113971.6357651352882385, 2220683829.8145909309387207, -370113971.6357656121253967, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357651352882385, -185056985.8178825676441193, -370113971.6357656121253967, 2624170957.4536480903625488, -E_100, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, E_100, -1034500000.0, -E_100, E_100, -740227943.2715302705764771, 370113971.6357653141021729, -1034500000.0, 1774727943.2715301513671875, -370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357653141021729, 185056985.8178827166557312};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
    }


    /**
    * Checks if the block builder and solver performs correctly the assemble of the extended system with constraints
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosExtendedDisplacementBlockBuilderAndSolverWithConstraints, KratosTrilinosApplicationMPITestSuite)
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
        ExtendedTestBuilderAndSolverDisplacement(r_model_part, r_comm, true);

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 22);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 22);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 204);

        // Values to check
        auto p_prop = r_model_part.pGetProperties(1);
        const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
        std::vector<int> row_indexes = {0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21};
        std::vector<int> column_indexes = {0, 1, 2, 4, 5, 9, 1, 2, 4, 5, 3, 1, 2, 4, 5, 8, 9, 12, 13, 1, 2, 4, 5, 8, 9, 11, 12, 13, 6, 7, 4, 5, 8, 9, 10, 1, 4, 5, 8, 9, 8, 10, 11, 12, 13, 14, 5, 10, 11, 12, 13, 4, 5, 10, 11, 12, 13, 16, 17, 4, 5, 10, 11, 12, 13, 15, 16, 17, 10, 14, 15, 16, 17, 18, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17, 20, 21, 12, 13, 14, 15, 16, 17, 19, 20, 21, 14, 18, 20, 17, 19, 16, 17, 18, 20, 21, 16, 17, 20, 21};
        std::vector<double> values = {740227943.2715302705764771, 1486220957.4536478519439697, -185056985.8178826570510864, 370113971.6357653141021729, -185056985.8178827166557312, -517250000.0, -185056985.8178827762603760, 1572984379.4520018100738525, -740227943.2715302705764771, 370113971.6357653141021729, 1.0, 370113971.6357653141021729, -740227943.2715302705764771, 1657021225.9261374473571777, -475379934.1969153881072998, -176565339.3830768167972565, -264848009.0746152102947235, -740227943.2715302705764771, 370113971.6357653141021729, -185056985.8178827166557312, 370113971.6357652544975281, -475379934.1969153881072998, 1457052651.9143548011779785, -264848009.0746152102947235, -397272013.6119228005409241, -689666666.6666666269302368, 370113971.6357653141021729, -185056985.8178827166557312, 1127028492.9089412689208984, 783913971.6357650756835938, -176565339.3830768167972565, -264848009.0746152102947235, 2245565339.3830766677856445, 264848009.0746151804924011, -1034500000.0, -517250000.0, -264848009.0746152102947235, -397272013.6119228005409241, 264848009.0746151804924011, 914522013.6119227409362793, -1034500000.0, 2434750982.5687417984008789, 365750982.5687416195869446, -365750982.5687417387962341, -365750982.5687416791915894, -1034500000.0, -689666666.6666666269302368, 365750982.5687416195869446, 1055417649.2354083061218262, -365750982.5687416791915894, -365750982.5687416195869446, -740227943.2715302705764771, 370113971.6357653141021729, -365750982.5687417387962341, -365750982.5687416791915894, 1846206869.1118023395538330, -374476960.7027890682220459, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -365750982.5687416791915894, -365750982.5687416195869446, -374476960.7027890682220459, 1770364954.2045071125030518, -1034500000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2809227943.2715301513671875, 370113971.6357650756835938, -740227943.2715302705764771, -370113971.6357651352882385, -1034500000.0, -1034500000.0, 370113971.6357650756835938, 1219556985.8178825378417969, -370113971.6357651352882385, -185056985.8178825676441193, -740227943.2715302705764771, 370113971.6357653141021729, -740227943.2715302705764771, -370113971.6357651352882385, 2220683829.8145909309387207, -370113971.6357656121253967, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357651352882385, -185056985.8178825676441193, -370113971.6357656121253967, 2624170957.4536480903625488, -E_100, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, E_100, -1034500000.0, -E_100, E_100, -740227943.2715302705764771, 370113971.6357653141021729, -1034500000.0, 1774727943.2715301513671875, -370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357653141021729, 185056985.8178827166557312};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 22);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 22);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 23);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
        column_indexes = {0, 1, 2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
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
    //     KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 19);
    //     KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 19);
    //     KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 161);

    //     // Values to check
    //     auto p_prop = r_model_part.pGetProperties(1);
    //     const double E_100 = p_prop->GetValue(YOUNG_MODULUS)/100.0;
    //     std::vector<int> row_indexes = {0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18};
    //     std::vector<int> column_indexes = {0, 1, 2, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 6, 1, 2, 3, 4, 5, 6, 9, 10, 1, 2, 3, 4, 5, 6, 8, 9, 10, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 5, 7, 8, 9, 10, 11, 4, 7, 8, 9, 10, 3, 4, 7, 8, 9, 10, 13, 14, 3, 4, 7, 8, 9, 10, 12, 13, 14, 7, 11, 12, 13, 14, 15, 10, 11, 12, 13, 14, 9, 10, 11, 12, 13, 14, 17, 18, 9, 10, 11, 12, 13, 14, 16, 17, 18, 11, 15, 17, 14, 16, 13, 14, 15, 17, 18, 13, 14, 17, 18};
    //     std::vector<double> values = {598856985.8178827762603760, 370113971.6357653141021729, -185056985.8178827166557312, 370113971.6357653141021729, 1572984379.4520018100738525, -555170957.4536480903625488, -740227943.2715302705764771, 370113971.6357653141021729, -185056985.8178827166557312, -555170957.4536479711532593, 1257477943.2715306282043457, 370113971.6357653141021729, -185056985.8178827166557312, -517250000.0, -740227943.2715302705764771, 370113971.6357653141021729, 1657021225.9261374473571777, -475379934.1969153881072998, -176565339.3830768167972565, -264848009.0746152102947235, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -475379934.1969153881072998, 1457052651.9143548011779785, -264848009.0746152102947235, -397272013.6119228005409241, -689666666.6666666269302368, 370113971.6357653141021729, -185056985.8178827166557312, -176565339.3830768167972565, -264848009.0746152102947235, 2245565339.3830766677856445, 264848009.0746151804924011, -1034500000.0, -517250000.0, -264848009.0746152102947235, -397272013.6119228005409241, 264848009.0746151804924011, 914522013.6119227409362793, -1034500000.0, 2434750982.5687417984008789, 365750982.5687416195869446, -365750982.5687417387962341, -365750982.5687416791915894, -1034500000.0, -689666666.6666666269302368, 365750982.5687416195869446, 1055417649.2354083061218262, -365750982.5687416791915894, -365750982.5687416195869446, -740227943.2715302705764771, 370113971.6357653141021729, -365750982.5687417387962341, -365750982.5687416791915894, 1846206869.1118023395538330, -374476960.7027890682220459, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -365750982.5687416791915894, -365750982.5687416195869446, -374476960.7027890682220459, 1770364954.2045071125030518, -1034500000.0, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, 2809227943.2715301513671875, 370113971.6357650756835938, -740227943.2715302705764771, -370113971.6357651352882385, -1034500000.0, -1034500000.0, 370113971.6357650756835938, 1219556985.8178825378417969, -370113971.6357651352882385, -185056985.8178825676441193, -740227943.2715302705764771, 370113971.6357653141021729, -740227943.2715302705764771, -370113971.6357651352882385, 2220683829.8145909309387207, -370113971.6357656121253967, -740227943.2715302705764771, 370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357651352882385, -185056985.8178825676441193, -370113971.6357656121253967, 2624170957.4536480903625488, -E_100, 370113971.6357653141021729, -185056985.8178827166557312, -1034500000.0, E_100, -1034500000.0, -E_100, E_100, -740227943.2715302705764771, 370113971.6357653141021729, -1034500000.0, 1774727943.2715301513671875, -370113971.6357653141021729, 370113971.6357652544975281, -185056985.8178827166557312, -370113971.6357653141021729, 185056985.8178827166557312};

    //     // Check assembly
    //     TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);
    // }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the extended system with constraints
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosStructuralDisplacementBlockBuilderAndSolver, KratosTrilinosApplicationMPITestSuite)
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
        const bool check = StructuralTestBuilderAndSolverDisplacement(r_model_part, r_comm);
        if (!check) {
            return void();
        }

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 54);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 54);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 1522);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 13, 14, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 12, 13, 14, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
        const double const0 = 111111.111111;
        const double const1 = 222222.222222;
        const double const2 = 83333.3333333;
        const double const3 = 138888.8888888;
        const double const4 = 55555.55555555;
        const double const5 = 444444.4444444;
        const double const6 = 227222.2222222;
        const double const7 = 454444.4444444;
        const double const8 = 116111.1111111;
        std::vector<double> values = {const6, const6, const1, const1, const1, const1, const6, const6, const1, const7, const4, const4, -const4, const7, const2, const4, const2, -const8, -const2, -const3, -const2, const2, const5, -const2, -const0, const2, const4, -const2, -const3, const1, const1, const1, const4, const5, -const4, const4, const4, -const2, const5, -const2, -const3, const2, -const0, const2, const2, -const0, -const2, const5, const2, -const3, -const2, const4, const4, -const4, const7, const4, -const8, const2, -const3, const2, const7, -const2, const4, -const2, -const2, const4, const2, -const3, -const2, const5, const2, -const0, -const4, const4, const4, const5, -const3, -const2, -const0, -const2, const4, const2, const5, const2, -const2, -const3, const2, const4, -const2, -const0, const2, const5, const6, const6, const1, const1, const1, const1, const6, const6, const1, const1, const1, const1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 54);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 54);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 72);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
        column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 0, 1, 2, 6, 7, 8, 9, 10, 11, 18, 19, 20, 24, 25, 26, 30, 31, 32};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

    /**
    * Checks if the block builder and solver performs correctly the assemble of the extended system with constraints
    */
    KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosStructuralDisplacementBlockBuilderAndSolverInverted, KratosTrilinosApplicationMPITestSuite)
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
        const bool check = StructuralTestBuilderAndSolverDisplacement(r_model_part, r_comm, true);
        if (!check) {
            return void();
        }

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
        KRATOS_EXPECT_EQ(rA.NumGlobalRows(), 54);
        KRATOS_EXPECT_EQ(rA.NumGlobalCols(), 54);
        KRATOS_EXPECT_EQ(rA.NumGlobalNonzeros(), 1522);

        // Values to check
        std::vector<int> row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 13, 14, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
        std::vector<int> column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 12, 13, 14, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 9, 15, 18, 21, 10, 11, 16, 17, 19, 20, 22, 23, 10, 11, 16, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
        const double const0 = 111111.111111;
        const double const1 = 222222.222222;
        const double const2 = 83333.3333333;
        const double const3 = 138888.8888888;
        const double const4 = 55555.55555555;
        const double const5 = 444444.4444444;
        const double const6 = 227222.2222222;
        const double const7 = 454444.4444444;
        const double const8 = 116111.1111111;
        std::vector<double> values = {const6, const6, const1, const1, const1, const1, const6, const6, const1, const7, const4, const4, -const4, const7, const2, const4, const2, -const8, -const2, -const3, -const2, const2, const5, -const2, -const0, const2, const4, -const2, -const3, const1, const1, const1, const4, const5, -const4, const4, const4, -const2, const5, -const2, -const3, const2, -const0, const2, const2, -const0, -const2, const5, const2, -const3, -const2, const4, const4, -const4, const7, const4, -const8, const2, -const3, const2, const7, -const2, const4, -const2, -const2, const4, const2, -const3, -const2, const5, const2, -const0, -const4, const4, const4, const5, -const3, -const2, -const0, -const2, const4, const2, const5, const2, -const2, -const3, const2, const4, -const2, -const0, const2, const5, const6, const6, const1, const1, const1, const1, const6, const6, const1, const1, const1, const1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly
        TrilinosCPPTestUtilities::CheckSparseMatrix(rA, row_indexes, column_indexes, values);

        // Now checking relation T matrix
        const auto& r_T = p_builder_and_solver->GetConstraintRelationMatrix();

        // // To create the solution of reference
        // DebugLHS(r_T);

        KRATOS_EXPECT_EQ(r_T.NumGlobalRows(), 54);
        KRATOS_EXPECT_EQ(r_T.NumGlobalCols(), 54);
        KRATOS_EXPECT_EQ(r_T.NumGlobalNonzeros(), 72);

        // Values to check
        row_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
        column_indexes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 0, 1, 2, 6, 7, 8, 9, 10, 11, 18, 19, 20, 24, 25, 26, 30, 31, 32};
        values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // Check assembly T matrix
        TrilinosCPPTestUtilities::CheckSparseMatrix(r_T, row_indexes, column_indexes, values);
    }

}  // namespace Kratos::Testing.