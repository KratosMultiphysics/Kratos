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
                pnode2->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
                auto pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
                pnode3->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

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
    }

    /**
    * @brief It generates a truss structure with an expected solution with an element with null contribution
    */
    static inline void BasicTestBuilderAndSolverDisplacementWithZeroContribution(ModelPart& rModelPart)
    {
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
        for (auto& node : rModelPart.Nodes()) {
            node.AddDof(DISPLACEMENT_X, REACTION_X);
            node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        }

        /// Initialize elements
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& elem : rModelPart.Elements()) {
            elem.Initialize(r_process_info);
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

        // Fix X in first node
        if (rModelPart.HasNode(1)) {
            auto pnode = rModelPart.pGetNode(1);
            pnode->Fix(DISPLACEMENT_X);
        }
    }

    /**
    * @brief It generates a truss structure with an expected solution
    */
    static inline void ExtendedTestBuilderAndSolverDisplacement(ModelPart& rModelPart, const bool WithConstraint = false)
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
        if (world_size == 1) {
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
        } else {
            KRATOS_ERROR << "Pending to do the proper partitions" << std::endl;
        }

        /// Add dof
        for (auto& node : rModelPart.Nodes()) {
            node.AddDof(DISPLACEMENT_X, REACTION_X);
            node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        }

        /// Initialize elements
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& elem : rModelPart.Elements()) {
            elem.Initialize(r_process_info);
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
        // Fix X in first node
        if (rModelPart.HasNode(1)) {
            auto pnode = rModelPart.pGetNode(1);
            pnode->Fix(DISPLACEMENT_X);
        }
        // Fix X and Y in fourth node
        if (rModelPart.HasNode(4)) {
            auto pnode = rModelPart.pGetNode(1);
            pnode->Fix(DISPLACEMENT_X);
            pnode->Fix(DISPLACEMENT_Y);
        }

        // // TODO: Fix this once MPC work!
        // if (WithConstraint) {
        //     rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode1, DISPLACEMENT_Y, *pnode2, DISPLACEMENT_Y, 1.0, 0.0);
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

    static void CheckAssembly(
        const TrilinosSparseSpaceType::MatrixType& rA,
        const std::vector<int>& rRowIndexes,
        const std::vector<int>& rColumnIndexes,
        const std::vector<double>& rValues,
        const double Tolerance = 1e-8
        )
    {
        int row, column;
        double value;
        for (std::size_t counter = 0; counter < rRowIndexes.size(); ++counter) {
            row = rRowIndexes[counter];
            column = rColumnIndexes[counter];
            value = rValues[counter];
            for (int i = 0; i < rA.NumMyRows(); i++) {
                int numEntries; // Number of non-zero entries
                double* vals;   // Row non-zero values
                int* cols;      // Column indices of row non-zero values
                rA.ExtractMyRowView(i, numEntries, vals, cols);
                const int row_gid = rA.RowMap().GID(i);
                if (row == row_gid) {
                    int j;
                    for (j = 0; j < numEntries; j++) {
                        const int col_gid = rA.ColMap().GID(cols[j]);
                        if (col_gid == column) {
                            KRATOS_CHECK_RELATIVE_NEAR(value, vals[j], Tolerance)
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }

    // static void DebugLHS(
    //     const TrilinosSparseSpaceType::MatrixType& rA,
    //     const DataCommunicator& rDataCommunicator
    //     )
    // {
    //     const int world_size = rDataCommunicator.Size();
    //     KRATOS_ERROR_IF_NOT(world_size == 1) << "Debug must be done with one MPI core" << std::endl;
    //     std::cout << "\n        KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), " << rA.NumGlobalRows() << ");\n";
    //     std::cout << "        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), " << rA.NumGlobalCols() << ");\n";
    //     std::cout << "        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), " << rA.NumGlobalNonzeros() << ");\n";

    //     std::vector<int> row_indexes;
    //     std::vector<int> column_indexes;
    //     std::vector<double> values;
    //     for (int i = 0; i < rA.NumMyRows(); i++) {
    //         int numEntries; // Number of non-zero entries
    //         double* vals;   // Row non-zero values
    //         int* cols;      // Column indices of row non-zero values
    //         rA.ExtractMyRowView(i, numEntries, vals, cols);
    //         const int row_gid = rA.RowMap().GID(i);
    //         int j;
    //         for (j = 0; j < numEntries; j++) {
    //             const int col_gid = rA.ColMap().GID(cols[j]);
    //             if (std::abs(vals[j]) > 0.99) {
    //                 row_indexes.push_back(row_gid);
    //                 column_indexes.push_back(col_gid);
    //                 values.push_back(vals[j]);
    //             }
    //         }
    //     }
    //     std::cout << "\n        // Values to check\n";
    //     std::cout << "        std::vector<int> row_indexes = {";
    //     for(std::size_t i = 0; i < row_indexes.size() - 1; ++i) {
    //         std::cout << row_indexes[i] << ", ";
    //     }
    //     std::cout << row_indexes[row_indexes.size() - 1] << "};";
    //     std::cout << "\n        std::vector<int> column_indexes = {";
    //     for(std::size_t i = 0; i < column_indexes.size() - 1; ++i) {
    //         std::cout << column_indexes[i] << ", ";
    //     }
    //     std::cout << column_indexes[column_indexes.size() - 1] << "};";
    //     std::cout << "\n        std::vector<double> values = {";
    //     for(std::size_t i = 0; i < values.size() - 1; ++i) {
    //         std::cout << std::fixed;
    //         std::cout << std::setprecision(16);
    //         std::cout << values[i] << ", ";
    //     }
    //     std::cout << std::fixed;
    //     std::cout << std::setprecision(16);
    //     std::cout << values[values.size() - 1] << "};" << std::endl;
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
        // DebugLHS(rA, r_comm);

        // The solution check
        KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 28);

        // Values to check
        std::vector<int> row_indexes = {0,2,2,4,4};
        std::vector<int> column_indexes = {0,2,4,2,4};
        std::vector<double> values = {2069000000.0000000000000000, 4138000000.0000000000000000, -2069000000.0000000000000000, -2069000000.0000000000000000, 2069000000.0000000000000000};

        // Check assembly
        CheckAssembly(rA, row_indexes, column_indexes, values);
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

        // Set MPI coomunicator
        ModelPartCommunicatorUtilities::SetMPICommunicator(r_model_part, r_comm);

        // Basic build
        BasicTestBuilderAndSolverDisplacementWithZeroContribution(r_model_part);

        // Compute communicaton plan and fill communicator meshes correctly
        ParallelFillCommunicator(r_model_part, r_comm).Execute();

        // Create the solvers and things required
        auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
        auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
        auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosBlockBuilderAndSolverType(epetra_comm, 15, p_solver) );

        const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

        // // To create the solution of reference
        // DebugLHS(rA, r_comm);

        // The solution check
        KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 6);
        KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 28);

        // Values to check
        std::vector<int> row_indexes = {0, 2};
        std::vector<int> column_indexes = {0, 2};
        std::vector<double> values = {2069000000.0000000000000000, 2069000000.0000000000000000};

        // Check assembly
        CheckAssembly(rA, row_indexes, column_indexes, values);
    }

    // /**
    // * Checks if the elimination builder and solver performs correctly the assemble of the system
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

    //     // Set MPI coomunicator
    //     ModelPartCommunicatorUtilities::SetMPICommunicator(r_model_part, r_comm);

    //     // Basic build
    //     BasicTestBuilderAndSolverDisplacement(r_model_part);

    //     // Compute communicaton plan and fill communicator meshes correctly
    //     ParallelFillCommunicator(r_model_part, r_comm).Execute();

    //     // Create the solvers and things required
    //     auto p_scheme = TrilinosSchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticSchemeType() );
    //     auto p_solver = TrilinosLinearSolverType::Pointer( new AmgclMPISolverType() );
    //     auto p_builder_and_solver = TrilinosBuilderAndSolverType::Pointer( new TrilinosResidualBasedEliminationBuilderAndSolverType(epetra_comm, 15, p_solver) );

    //     const auto& rA = BuildSystem(r_model_part, p_scheme, p_builder_and_solver);

    //     // // To create the solution of reference
    //     // DebugLHS(rA, r_comm);

    //     // The solution check
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalRows(), 2);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalCols(), 2);
    //     KRATOS_CHECK_EQUAL(rA.NumGlobalNonzeros(), 4);

    //     // Values to check
    //     std::vector<int> row_indexes = {0, 0, 1, 1};
    //     std::vector<int> column_indexes = {0, 1, 0, 1};
    //     std::vector<double> values = {4138000000.0000000000000000, -2069000000.0000000000000000, -2069000000.0000000000000000, 2069000000.0000000000000000};

    //     // Check assembly
    //     CheckAssembly(rA, row_indexes, column_indexes, values);
    // }

}  // namespace Kratos::Testing.