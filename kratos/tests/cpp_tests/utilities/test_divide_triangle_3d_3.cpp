//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//                   Kratos default license: kratos/license.txt
//
//   Main authors:   Pablo Becker
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/expect.h"
#include "utilities/divide_triangle_3d_3.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle3D3, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, -1.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the triangle splitting utility
    DivideTriangle3D3<Node> triangle_splitter(r_geometry, r_elemental_distances);

    // Call the divide geometry method
    triangle_splitter.GenerateDivision();


    // Call the intersection generation method
    triangle_splitter.GenerateIntersectionsSkin();

    // Call the positive exterior faces generation method
    std::vector < unsigned int > pos_ext_faces_parent_ids;
    std::vector < DivideTriangle2D3<Node>::IndexedPointGeometryPointerType > pos_ext_faces;
    triangle_splitter.GenerateExteriorFaces(
        pos_ext_faces,
        pos_ext_faces_parent_ids,
        triangle_splitter.GetPositiveSubdivisions());

    // Call the negative exterior faces generation method
    std::vector < unsigned int > neg_ext_faces_parent_ids;
    std::vector < DivideTriangle2D3<Node>::IndexedPointGeometryPointerType > neg_ext_faces;
    triangle_splitter.GenerateExteriorFaces(
        neg_ext_faces,
        neg_ext_faces_parent_ids,
        triangle_splitter.GetNegativeSubdivisions());

    const double tolerance = 1e-10;

    // Check general splitting values
    KRATOS_EXPECT_TRUE(triangle_splitter.mIsSplit);
    KRATOS_EXPECT_EQ(triangle_splitter.mDivisionsNumber, 3);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdgesNumber, 2);

    // Check split edges
    std::array<int,6> expected_split_edges;
    expected_split_edges[0]= 0;
    expected_split_edges[1]= 1;
    expected_split_edges[2]= 2;
    expected_split_edges[3]= -1;
    expected_split_edges[4]= 4;
    expected_split_edges[5]= 5;
    for (unsigned int i = 0; i < 6; ++i) {
        KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[i],  expected_split_edges[i]);
    }

    // Check subdivisions
    const auto &r_positive_subdivision_0 = *(triangle_splitter.GetPositiveSubdivisions()[0]);

    Vector positive_subdivision_0_node_0_coordinates(3);
    positive_subdivision_0_node_0_coordinates[0] = 0.0;
    positive_subdivision_0_node_0_coordinates[1] = 0.5;
    positive_subdivision_0_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_positive_subdivision_0[0],positive_subdivision_0_node_0_coordinates);

    Vector positive_subdivision_0_node_1_coordinates(3);
    positive_subdivision_0_node_1_coordinates[0] = 0.5;
    positive_subdivision_0_node_1_coordinates[1] = 0.5;
    positive_subdivision_0_node_1_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_positive_subdivision_0[1],positive_subdivision_0_node_1_coordinates);

    Vector positive_subdivision_0_node_2_coordinates(3);
    positive_subdivision_0_node_2_coordinates[0] = 0.0;
    positive_subdivision_0_node_2_coordinates[1] = 1.0;
    positive_subdivision_0_node_2_coordinates[2] = -1.0;
    KRATOS_EXPECT_VECTOR_EQ(r_positive_subdivision_0[2],positive_subdivision_0_node_2_coordinates);

    KRATOS_EXPECT_NEAR(r_positive_subdivision_0.Area(), 0.17677669529,tolerance);

    const auto &r_negative_subdivision_0 = *(triangle_splitter.GetNegativeSubdivisions()[0]);

    Vector negative_subdivision_0_node_0_coordinates(3);
    negative_subdivision_0_node_0_coordinates[0] = 0.0;
    negative_subdivision_0_node_0_coordinates[1] = 0.5;
    negative_subdivision_0_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_subdivision_0[0],negative_subdivision_0_node_0_coordinates);

    Vector negative_subdivision_0_node_1_coordinates(3);
    negative_subdivision_0_node_1_coordinates[0] = 1.0;
    negative_subdivision_0_node_1_coordinates[1] = 0.0;
    negative_subdivision_0_node_1_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_subdivision_0[1],negative_subdivision_0_node_1_coordinates);

    Vector negative_subdivision_0_node_2_coordinates(3);
    negative_subdivision_0_node_2_coordinates[0] = 0.5;
    negative_subdivision_0_node_2_coordinates[1] = 0.5;
    negative_subdivision_0_node_2_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_subdivision_0[2],negative_subdivision_0_node_2_coordinates);

    KRATOS_EXPECT_NEAR(r_negative_subdivision_0.Area(), 0.17677669529, tolerance);

    const auto &r_negative_subdivision_1 = *(triangle_splitter.GetNegativeSubdivisions()[1]);

    Vector negative_subdivision_1_node_0_coordinates(3);
    negative_subdivision_1_node_0_coordinates[0] = 0.0;
    negative_subdivision_1_node_0_coordinates[1] = 0.5;
    negative_subdivision_1_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_subdivision_1[0],negative_subdivision_1_node_0_coordinates);

    Vector negative_subdivision_1_node_1_coordinates(3);
    negative_subdivision_1_node_1_coordinates[0] = 0.0;
    negative_subdivision_1_node_1_coordinates[1] = 0.0;
    negative_subdivision_1_node_1_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_subdivision_1[1],negative_subdivision_1_node_1_coordinates);

    Vector negative_subdivision_1_node_2_coordinates(3);
    negative_subdivision_1_node_2_coordinates[0] = 1.0;
    negative_subdivision_1_node_2_coordinates[1] = 0.0;
    negative_subdivision_1_node_2_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_subdivision_1[2],negative_subdivision_1_node_2_coordinates);

    KRATOS_EXPECT_NEAR(r_negative_subdivision_1.Area(), 0.35355339059, tolerance);

    // Check interfaces
    const auto &r_positive_interface_0 = *(triangle_splitter.GetPositiveInterfaces()[0]);

    Vector positive_interface_0_node_0_coordinates(3);
    positive_interface_0_node_0_coordinates[0] = 0.0;
    positive_interface_0_node_0_coordinates[1] = 0.5;
    positive_interface_0_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_positive_interface_0[0],positive_interface_0_node_0_coordinates);

    Vector positive_interface_0_node_1_coordinates(3);
    positive_interface_0_node_1_coordinates[0] = 0.5;
    positive_interface_0_node_1_coordinates[1] = 0.5;
    positive_interface_0_node_1_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_positive_interface_0[1],positive_interface_0_node_1_coordinates);

    const auto &r_negative_interface_0 = *(triangle_splitter.GetNegativeInterfaces()[0]);

    Vector negative_interface_0_node_0_coordinates(3);
    negative_interface_0_node_0_coordinates[0] = 0.5;
    negative_interface_0_node_0_coordinates[1] = 0.5;
    negative_interface_0_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_interface_0[0],negative_interface_0_node_0_coordinates);

    Vector negative_interface_0_node_1_coordinates(3);
    negative_interface_0_node_1_coordinates[0] = 0.0;
    negative_interface_0_node_1_coordinates[1] = 0.5;
    negative_interface_0_node_1_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ(r_negative_interface_0[1],negative_interface_0_node_1_coordinates);

    KRATOS_EXPECT_EQ(triangle_splitter.GetPositiveInterfacesParentIds()[0], 0);
    KRATOS_EXPECT_EQ(triangle_splitter.GetNegativeInterfacesParentIds()[0], 0);

    // Check exterior faces
    KRATOS_EXPECT_EQ(pos_ext_faces.size(), 2);
    KRATOS_EXPECT_EQ(neg_ext_faces.size(), 3);

    KRATOS_EXPECT_EQ(pos_ext_faces_parent_ids[0], 0);
    KRATOS_EXPECT_EQ(pos_ext_faces_parent_ids[1], 0);

    KRATOS_EXPECT_EQ(neg_ext_faces_parent_ids[0], 0);
    KRATOS_EXPECT_EQ(neg_ext_faces_parent_ids[1], 1);
    KRATOS_EXPECT_EQ(neg_ext_faces_parent_ids[2], 1);

    Vector pos_ext_faces_0_node_0_coordinates(3);
    pos_ext_faces_0_node_0_coordinates[0] = 0.5;
    pos_ext_faces_0_node_0_coordinates[1] = 0.5;
    pos_ext_faces_0_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ((*pos_ext_faces[0])[0],pos_ext_faces_0_node_0_coordinates);

    Vector pos_ext_faces_0_node_1_coordinates(3);
    pos_ext_faces_0_node_1_coordinates[0] = 0.0;
    pos_ext_faces_0_node_1_coordinates[1] = 1.0;
    pos_ext_faces_0_node_1_coordinates[2] = -1.0;
    KRATOS_EXPECT_VECTOR_EQ((*pos_ext_faces[0])[1],pos_ext_faces_0_node_1_coordinates);

    Vector pos_ext_faces_1_node_0_coordinates(3);
    pos_ext_faces_1_node_0_coordinates[0] = 0.0;
    pos_ext_faces_1_node_0_coordinates[1] = 1.0;
    pos_ext_faces_1_node_0_coordinates[2] = -1.0;
    KRATOS_EXPECT_VECTOR_EQ((*pos_ext_faces[1])[0],pos_ext_faces_1_node_0_coordinates);

    Vector pos_ext_faces_1_node_1_coordinates(3);
    pos_ext_faces_1_node_1_coordinates[0] = 0.0;
    pos_ext_faces_1_node_1_coordinates[1] = 0.5;
    pos_ext_faces_1_node_1_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ((*pos_ext_faces[1])[1],pos_ext_faces_1_node_1_coordinates);

    Vector neg_ext_faces_0_node_0_coordinates(3);
    neg_ext_faces_0_node_0_coordinates[0] = 1.0;
    neg_ext_faces_0_node_0_coordinates[1] = 0.0;
    neg_ext_faces_0_node_0_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ((*neg_ext_faces[0])[0],neg_ext_faces_0_node_0_coordinates);

    Vector neg_ext_faces_0_node_1_coordinates(3);
    neg_ext_faces_0_node_1_coordinates[0] = 0.5;
    neg_ext_faces_0_node_1_coordinates[1] = 0.5;
    neg_ext_faces_0_node_1_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ((*neg_ext_faces[0])[1],neg_ext_faces_0_node_1_coordinates);

    Vector neg_ext_faces_1_node_0_coordinates(3);
    neg_ext_faces_1_node_0_coordinates[0] = 0.0;
    neg_ext_faces_1_node_0_coordinates[1] = 0.5;
    neg_ext_faces_1_node_0_coordinates[2] = -0.5;
    KRATOS_EXPECT_VECTOR_EQ((*neg_ext_faces[1])[0],neg_ext_faces_1_node_0_coordinates);

    Vector neg_ext_faces_1_node_1_coordinates(3);
    neg_ext_faces_1_node_1_coordinates[0] = 0.0;
    neg_ext_faces_1_node_1_coordinates[1] = 0.0;
    neg_ext_faces_1_node_1_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ((*neg_ext_faces[1])[1],neg_ext_faces_1_node_1_coordinates);

    Vector neg_ext_faces_2_node_0_coordinates(3);
    neg_ext_faces_2_node_0_coordinates[0] = 0.0;
    neg_ext_faces_2_node_0_coordinates[1] = 0.0;
    neg_ext_faces_2_node_0_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ((*neg_ext_faces[2])[0],neg_ext_faces_2_node_0_coordinates);

    Vector neg_ext_faces_2_node_1_coordinates(3);
    neg_ext_faces_2_node_1_coordinates[0] = 1.0;
    neg_ext_faces_2_node_1_coordinates[1] = 0.0;
    neg_ext_faces_2_node_1_coordinates[2] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ((*neg_ext_faces[2])[1],neg_ext_faces_2_node_1_coordinates);

}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle3D3ZeroNode, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.3,  0.0, 0.6);
    base_model_part.CreateNewNode(2, 0.3,  0.0, 0.9);
    base_model_part.CreateNewNode(3, 0.3, -0.3, 0.3);
    Properties::Pointer p_properties(new Properties(0));
    auto p_elem = base_model_part.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties);

    Vector nodal_distances = ZeroVector(3);
    // cut at z = 0.6, across node 0
    nodal_distances[0] = 0.0;
    nodal_distances[1] = 0.3;
    nodal_distances[2] = -0.3;

    auto divider = DivideTriangle3D3<Node>(p_elem->GetGeometry(), nodal_distances);
    divider.GenerateDivision();

    // Should split a single edge and return two triangles
    KRATOS_EXPECT_EQ(divider.GetPositiveSubdivisions().size(), 1);
    KRATOS_EXPECT_EQ(divider.GetNegativeSubdivisions().size(), 1);
}

}  // namespace Kratos::Testing.
