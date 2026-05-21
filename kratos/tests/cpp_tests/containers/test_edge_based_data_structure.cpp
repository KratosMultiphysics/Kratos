//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/edge_based_data_structure.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
//#include "includes/gid_io.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "testing/testing.h"

namespace Kratos::Testing {

namespace
{
    // Auxiliary function to check a vector of indices
    // TODO: Remove as soon as we update the testing
    void CheckIndicesVector(
        const std::vector<std::size_t>& rVector,
        const std::vector<std::size_t>& rRefVector)
    {
        KRATOS_CHECK(rVector.size() == rRefVector.size());
        for (std::size_t i = 0; i < rVector.size(); ++i) {
            KRATOS_CHECK(rVector[i] == rRefVector[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(EdgeBasedDataStructure2D, KratosCoreFastSuite)
{
    // Set-up a simplicial mesh to calculate its edge data structure
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");
    auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions" : 2,
        "element_name" : "Element2D3N",
        "condition_name" : "LineCondition",
        "create_skin_sub_model_part" : true
    })");
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

    //GidIO<> gid_io_convection("/home/rzorrilla/Desktop/EdgeBasedDataStructure2D", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    //gid_io_convection.InitializeMesh(0);
    //gid_io_convection.WriteMesh(r_model_part.GetMesh());
    //gid_io_convection.FinalizeMesh();

    // Calculate nodal neighbours
    // TODO: Temporary solution until we skip the neighbours calculation
    FindGlobalNodalNeighboursProcess nodal_neighs_process(r_model_part);
    nodal_neighs_process.Execute();

    // Create the edge-based data structure
    EdgeBasedDataStructure<2> edge_based_data_structure;
    edge_based_data_structure.CalculateEdgeDataStructure(r_model_part);

    // Check edge 1-2
    auto &r_edge_data_12 = edge_based_data_structure.GetEdgeData(1, 2);
    KRATOS_CHECK(r_edge_data_12.IsBoundary())
    KRATOS_CHECK_NEAR(r_edge_data_12.GetLength(), 0.5, 1.0e-8);
    KRATOS_CHECK_NEAR(r_edge_data_12.GetOffDiagonalConsistentMass(), 0.01041666, 1.0e-8);
    KRATOS_CHECK_NEAR(r_edge_data_12.GetOffDiagonalLaplacian(), -0.5, 1.0e-8);
    const std::vector<double> off_diag_conv_ref_12({-0.083333333, 0.083333333});
    const std::vector<double> off_diag_conv_trans_ref_12({0.0, -0.083333333});
    const std::vector<double> off_diag_conv_boundary_ref_12({-0.125, 0.0});
    KRATOS_CHECK_VECTOR_NEAR(r_edge_data_12.GetOffDiagonalConvective(), off_diag_conv_ref_12, 1.0e-8)
    KRATOS_CHECK_VECTOR_NEAR(r_edge_data_12.GetOffDiagonalConvectiveTranspose(), off_diag_conv_trans_ref_12, 1.0e-8)
    KRATOS_CHECK_VECTOR_NEAR(r_edge_data_12.GetOffDiagonalConvectiveBoundary(), off_diag_conv_boundary_ref_12, 1.0e-8)

    // Check edge 1-5
    auto &r_edge_data_15 = edge_based_data_structure.GetEdgeData(1, 5);
    KRATOS_CHECK_IS_FALSE(r_edge_data_15.IsBoundary())
    KRATOS_CHECK_NEAR(r_edge_data_15.GetLength(), 0.5*std::sqrt(2), 1.0e-8);
    KRATOS_CHECK_NEAR(r_edge_data_15.GetOffDiagonalConsistentMass(), 0.020833333, 1.0e-8);
    KRATOS_CHECK_NEAR(r_edge_data_15.GetOffDiagonalLaplacian(), 0.0, 1.0e-8);
    const std::vector<double> off_diag_conv_ref_15({0.083333333, 0.083333333});
    const std::vector<double> off_diag_conv_trans_ref_15({-0.083333333, -0.083333333});
    KRATOS_CHECK_VECTOR_NEAR(r_edge_data_15.GetOffDiagonalConvective(), off_diag_conv_ref_15, 1.0e-8)
    KRATOS_CHECK_VECTOR_NEAR(r_edge_data_15.GetOffDiagonalConvectiveTranspose(), off_diag_conv_trans_ref_15, 1.0e-8)

    // Check diagonal components
    const std::vector<double> bd_mass_mat_diag_1_ref({-0.125,-0.125});
    const std::vector<double> mass_mat_diag_ref({0.04166666, 0.0625, 0.020833333, 0.0625, 0.125, 0.0625, 0.020833333, 0.0625, 0.04166666});
    const std::vector<double> lumped_mass_mat_diag_ref({0.083333333, 0.125, 0.04166666, 0.125, 0.25, 0.125, 0.04166666, 0.125, 0.083333333});
    KRATOS_CHECK_VECTOR_NEAR(edge_based_data_structure.GetMassMatrixDiagonal(), mass_mat_diag_ref, 1.0e-8);
    KRATOS_CHECK_VECTOR_NEAR(edge_based_data_structure.GetLumpedMassMatrixDiagonal(), lumped_mass_mat_diag_ref, 1.0e-8);
    KRATOS_CHECK_VECTOR_NEAR(edge_based_data_structure.GetBoundaryMassMatrixDiagonal(1), bd_mass_mat_diag_1_ref, 1.0e-8);

    // Check sparse data
    const std::vector<std::size_t> row_indices_ref({0,0,3,6,7,10,13,14,15,16});
    const std::vector<std::size_t> col_indices_ref({2,4,5,3,5,6,6,5,7,8,6,8,9,9,8,9});
    CheckIndicesVector(edge_based_data_structure.GetRowIndices(), row_indices_ref);
    CheckIndicesVector(edge_based_data_structure.GetColIndices(), col_indices_ref);
    KRATOS_CHECK_EQUAL(edge_based_data_structure.NumberOfEdges(), 16);
}

} // namespace Kratos::Testing.
