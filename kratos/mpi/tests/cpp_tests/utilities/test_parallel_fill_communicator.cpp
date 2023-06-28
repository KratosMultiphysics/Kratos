//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/line_2d_2.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_bar_element.h"

namespace Kratos::Testing 
{

/**
* @brief It generates a truss structure with an expected solution
*/
static inline void GenerateModelPartEntinties(
    ModelPart& rModelPart,
    const DataCommunicator& rDataCommunicator
    )
{
    // Set MPI coomunicator
    ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, rDataCommunicator);

    // Add variables
    rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Create properties
    auto p_prop = rModelPart.CreateNewProperties(1, 0);

    // MPI data
    const int rank = rDataCommunicator.Rank();
    const int world_size = rDataCommunicator.Size();

    // Initially everything in one partition
    if (world_size == 1) {
        auto pnode1 = rModelPart.CreateNewNode(1, 10.0, -5.0, 0.0);
        auto pnode2 = rModelPart.CreateNewNode(2, 8.0, -4.0, 0.0);
        auto pnode3 = rModelPart.CreateNewNode(3, 6.0, -3.0, 0.0);
        auto pnode4 = rModelPart.CreateNewNode(4, 10.0, 0.0, 0.0);
        auto pnode5 = rModelPart.CreateNewNode(5, 8.0, 0.0, 0.0);
        auto pnode6 = rModelPart.CreateNewNode(6, 6.0, 0.0, 0.0);
        auto pnode7 = rModelPart.CreateNewNode(7, 4.0, -2.0, 0.0);
        auto pnode8 = rModelPart.CreateNewNode(8, 4.0, 0.0, 0.0);
        auto pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
        auto pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);
        auto pnode11 = rModelPart.CreateNewNode(11, 0.0, 0.0, 0.0);

        /// Add PARTITION_INDEX
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
        }

        auto pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode11, pnode10})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
        auto pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode10, pnode8})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
        auto pgeom3 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode8, pnode6})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 3, pgeom3, p_prop));
        auto pgeom4 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode6, pnode5})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 4, pgeom4, p_prop));
        auto pgeom5 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode5, pnode4})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 5, pgeom5, p_prop));
        auto pgeom6 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode4, pnode1})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 6, pgeom6, p_prop));
        auto pgeom7 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 7, pgeom7, p_prop));
        auto pgeom8 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 8, pgeom8, p_prop));
        auto pgeom9 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode3, pnode7})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 9, pgeom9, p_prop));
        auto pgeom10 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode7, pnode9})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 10, pgeom10, p_prop));
        auto pgeom11 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode9, pnode11})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 11, pgeom11, p_prop));
        auto pgeom12 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode10, pnode9})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 12, pgeom12, p_prop));
        auto pgeom13 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode9, pnode8})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 13, pgeom13, p_prop));
        auto pgeom14 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode8, pnode7})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 14, pgeom14, p_prop));
        auto pgeom15 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode7, pnode6})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 15, pgeom15, p_prop));
        auto pgeom16 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode6, pnode3})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 16, pgeom16, p_prop));
        auto pgeom17 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode3, pnode5})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 17, pgeom17, p_prop));
        auto pgeom18 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode5, pnode2})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 18, pgeom18, p_prop));
        auto pgeom19 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode4})});
        rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 19, pgeom19, p_prop));
    } else { // if (world_size == 1) { // TODO: Do more than one partition
        if (rank == 0) {
            auto pnode1 = rModelPart.CreateNewNode(1, 10.0, -5.0, 0.0);
            auto pnode2 = rModelPart.CreateNewNode(2, 8.0, -4.0, 0.0);
            auto pnode3 = rModelPart.CreateNewNode(3, 6.0, -3.0, 0.0);
            auto pnode4 = rModelPart.CreateNewNode(4, 10.0, 0.0, 0.0);
            auto pnode5 = rModelPart.CreateNewNode(5, 8.0, 0.0, 0.0);
            auto pnode6 = rModelPart.CreateNewNode(6, 6.0, 0.0, 0.0);
            auto pnode7 = rModelPart.CreateNewNode(7, 4.0, -2.0, 0.0);
            auto pnode8 = rModelPart.CreateNewNode(8, 4.0, 0.0, 0.0);
            auto pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
            auto pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }
            pnode9->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            pnode10->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

            auto pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode10, pnode8})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 2, pgeom2, p_prop));
            auto pgeom3 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode8, pnode6})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 3, pgeom3, p_prop));
            auto pgeom4 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode6, pnode5})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 4, pgeom4, p_prop));
            auto pgeom5 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode5, pnode4})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 5, pgeom5, p_prop));
            auto pgeom6 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode4, pnode1})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 6, pgeom6, p_prop));
            auto pgeom7 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 7, pgeom7, p_prop));
            auto pgeom8 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 8, pgeom8, p_prop));
            auto pgeom9 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode3, pnode7})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 9, pgeom9, p_prop));
            auto pgeom10 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode7, pnode9})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 10, pgeom10, p_prop));
            auto pgeom12 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode10, pnode9})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 12, pgeom12, p_prop));
            auto pgeom13 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode9, pnode8})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 13, pgeom13, p_prop));
            auto pgeom14 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode8, pnode7})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 14, pgeom14, p_prop));
            auto pgeom15 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode7, pnode6})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 15, pgeom15, p_prop));
            auto pgeom16 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode6, pnode3})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 16, pgeom16, p_prop));
            auto pgeom17 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode3, pnode5})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 17, pgeom17, p_prop));
            auto pgeom18 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode5, pnode2})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 18, pgeom18, p_prop));
            auto pgeom19 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode4})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 19, pgeom19, p_prop));
        } else if (rank == 1) {
            auto pnode9 = rModelPart.CreateNewNode(9, 2.0, -1.0, 0.0);
            auto pnode10 = rModelPart.CreateNewNode(10, 2.0, 0.0, 0.0);
            auto pnode11 = rModelPart.CreateNewNode(11, 0.0, 0.0, 0.0);

            /// Add PARTITION_INDEX
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = rank;
            }

            auto pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode11, pnode10})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 1, pgeom1, p_prop));
            auto pgeom11 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode9, pnode11})});
            rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>( 11, pgeom11, p_prop));
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(ParallelFillCommunicatorExecute, KratosMPICoreFastSuite)
{
    // The model part
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    
    // The data communicator
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();
    
    // MPI data
    const int rank =  r_data_communicator.Rank();
    const int world_size = r_data_communicator.Size();

    // Fill the model part
    GenerateModelPartEntinties(r_model_part, r_data_communicator);

    // Compute communicaton plan and fill communicator meshes correctly
    auto filler = ParallelFillCommunicator(r_model_part, r_data_communicator);
    filler.Execute();

    // // Debug info
    // filler.PrintModelPartDebugInfo(r_model_part);

    // Check that the communicator is correctly filled
    const auto& r_neighbours_indices = r_model_part.GetCommunicator().NeighbourIndices();
    if (world_size == 1) {
        KRATOS_CHECK_EQUAL(r_neighbours_indices.size(), 0);
        KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 11);
        KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
        KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 0);
    } else if (world_size == 2) {
        KRATOS_CHECK_EQUAL(r_neighbours_indices.size(), 1);
        if (rank == 0) {
            KRATOS_CHECK_EQUAL(r_neighbours_indices[0], 1);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 8);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 2);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        } else if (rank == 1) {
            KRATOS_CHECK_EQUAL(r_neighbours_indices[0], 0);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 3);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        }
    } else {
        if (rank == 0) {
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 8);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 2);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        } else if (rank == 1) {
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 3);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 2);
        } else { // The rest of the ranks
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().LocalMesh().NumberOfNodes(), 0);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GhostMesh().NumberOfNodes(), 0);
            KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().InterfaceMesh().NumberOfNodes(), 0);
        }
    }
}

} // namespace Kratos::Testing