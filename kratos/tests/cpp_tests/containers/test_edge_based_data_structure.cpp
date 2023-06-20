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
#include "includes/gid_io.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "testing/testing.h"

namespace Kratos::Testing {

void PrintVector(const std::vector<double>& rVect)
{
    std::cout << std::endl;
    for (const double& val : rVect) {
        std::cout << val << std::endl;
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

    GidIO<> gid_io_convection("/home/rzorrilla/Desktop/EdgeBasedDataStructure2D", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    gid_io_convection.InitializeMesh(0);
    gid_io_convection.WriteMesh(r_model_part.GetMesh());
    gid_io_convection.FinalizeMesh();

    // Calculate nodal neighbours
    // TODO: Temporary solution until we skip the neighbours calculation
    FindGlobalNodalNeighboursProcess nodal_neighs_process(r_model_part);
    nodal_neighs_process.Execute();

    // Create the edge-based data structure
    EdgeBasedDataStructure<2> edge_based_data_structure;
    edge_based_data_structure.CalculateEdgeDataStructure(r_model_part);

    auto &r_edge_data_12 = edge_based_data_structure.GetEdgeData(1, 2);
    KRATOS_WATCH(r_edge_data_12.GetOffDiagonalConsistentMass())
    KRATOS_WATCH(r_edge_data_12.GetOffDiagonalLaplacian())
    KRATOS_WATCH(r_edge_data_12.GetOffDiagonalConvective())
    KRATOS_WATCH(r_edge_data_12.GetOffDiagonalConvectiveTranspose())

    auto &r_edge_data_15 = edge_based_data_structure.GetEdgeData(1, 5);
    KRATOS_WATCH(r_edge_data_15.GetOffDiagonalConsistentMass())
    KRATOS_WATCH(r_edge_data_15.GetOffDiagonalLaplacian())
    KRATOS_WATCH(r_edge_data_15.GetOffDiagonalConvective())
    KRATOS_WATCH(r_edge_data_15.GetOffDiagonalConvectiveTranspose())

    const auto& r_mass_mat_diag = edge_based_data_structure.GetMassMatrixDiagonal();
    const auto& r_lumped_mass_mat_diag = edge_based_data_structure.GetLumpedMassMatrixDiagonal();
    PrintVector(r_mass_mat_diag);
    PrintVector(r_lumped_mass_mat_diag);

    KRATOS_CHECK_EQUAL(edge_based_data_structure.NumberOfEdges(), 16);

}

} // namespace Kratos::Testing.
