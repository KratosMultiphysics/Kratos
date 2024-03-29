//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/expect.h"
#include "utilities/divide_triangle_2d_3.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3Horizontal, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

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
    DivideTriangle2D3<Node> triangle_splitter(r_geometry, r_elemental_distances);

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
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[0],  0);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[1],  1);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[2],  2);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[3], -1);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[4],  4);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[5],  5);

    // Check subdivisions
    const auto &r_positive_subdivision_0 = *(triangle_splitter.GetPositiveSubdivisions()[0]);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[1].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[2].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[2].Y(), 1.0, tolerance);

    const auto &r_negative_subdivision_0 = *(triangle_splitter.GetNegativeSubdivisions()[0]);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[1].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[1].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[2].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[2].Y(), 0.5, tolerance);

    const auto &r_negative_subdivision_1 = *(triangle_splitter.GetNegativeSubdivisions()[1]);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[1].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[2].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[2].Y(), 0.0, tolerance);

    // Check interfaces
    const auto &r_positive_interface_0 = *(triangle_splitter.GetPositiveInterfaces()[0]);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[1].Y(), 0.5, tolerance);

    const auto &r_negative_interface_0 = *(triangle_splitter.GetNegativeInterfaces()[0]);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[1].Y(), 0.5, tolerance);

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

    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[1].Y(), 1.0, tolerance);

    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[0].Y(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[1].Y(), 0.5, tolerance);

    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[0].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[1].Y(), 0.5, tolerance);

    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[1].Y(), 0.0, tolerance);

    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[1].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[1].Y(), 0.0, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3Vertical, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) =  1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = -1.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the triangle splitting utility
    DivideTriangle2D3<Node> triangle_splitter(r_geometry, r_elemental_distances);

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
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[0],  0);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[1],  1);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[2],  2);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[3],  3);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[4],  4);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdges[5], -1);

    // Check subdivisions
    const auto &r_positive_subdivision_0 = *(triangle_splitter.GetPositiveSubdivisions()[0]);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[1].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[2].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_subdivision_0[2].Y(), 0.0, tolerance);

    const auto &r_negative_subdivision_0 = *(triangle_splitter.GetNegativeSubdivisions()[0]);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[1].Y(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[2].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_0[2].Y(), 0.0, tolerance);

    const auto &r_negative_subdivision_1 = *(triangle_splitter.GetNegativeSubdivisions()[1]);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[1].Y(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[2].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_subdivision_1[2].Y(), 0.0, tolerance);

    // Check interfaces
    const auto &r_positive_interface_0 = *(triangle_splitter.GetPositiveInterfaces()[0]);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_positive_interface_0[1].Y(), 0.0, tolerance);
    const auto &r_negative_interface_0 = *(triangle_splitter.GetNegativeInterfaces()[0]);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR(r_negative_interface_0[1].Y(), 0.5, tolerance);

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

    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[0].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[0])[1].Y(), 0.5, tolerance);

    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[1].X(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR((*pos_ext_faces[1])[1].Y(), 0.0, tolerance);

    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[0].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[0].Y(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[0])[1].Y(), 1.0, tolerance);

    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[0].Y(), 1.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[1])[1].Y(), 0.0, tolerance);

    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[0].X(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[0].Y(), 0.0, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[1].X(), 0.5, tolerance);
    KRATOS_EXPECT_NEAR((*neg_ext_faces[2])[1].Y(), 0.0, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3NoDivision, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < r_geometry.PointsNumber(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the triangle splitting utility
    DivideTriangle2D3<Node> triangle_splitter(r_geometry, r_elemental_distances);

    // Call the divide geometry method
    triangle_splitter.GenerateDivision();

    // Check general splitting values
    KRATOS_EXPECT_FALSE(triangle_splitter.mIsSplit);
    KRATOS_EXPECT_EQ(triangle_splitter.mDivisionsNumber, 1);
    KRATOS_EXPECT_EQ(triangle_splitter.mSplitEdgesNumber, 0);

}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3ZeroNode, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.6, 0.0);
    base_model_part.CreateNewNode(2, 0.0, 0.9, 0.0);
    base_model_part.CreateNewNode(3,-0.3, 0.3, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    auto p_elem = base_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

    Vector nodal_distances = ZeroVector(3);
    // cut at y = 0.6, across node 0
    nodal_distances[0] = 0.0;
    nodal_distances[1] = 0.3;
    nodal_distances[2] = -0.3;

    auto divider = DivideTriangle2D3<Node>(p_elem->GetGeometry(), nodal_distances);
    divider.GenerateDivision();

    // Should split a single edge and return two triangles
    KRATOS_EXPECT_EQ(divider.GetPositiveSubdivisions().size(), 1);
    KRATOS_EXPECT_EQ(divider.GetNegativeSubdivisions().size(), 1);
}


KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3TwoZeroNodes, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.6, 0.0);
    base_model_part.CreateNewNode(2, 0.0, 0.9, 0.0);
    base_model_part.CreateNewNode(3,-0.3, 0.3, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    auto p_elem = base_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

    Vector nodal_distances = ZeroVector(3);
    // cut at y = 0.6, across node 0
    nodal_distances[0] = 0.0;
    nodal_distances[1] = 0.0;
    nodal_distances[2] = -0.3;

    auto divider_one_negative = DivideTriangle2D3<Node>(p_elem->GetGeometry(), nodal_distances);
    divider_one_negative.GenerateDivision();

    KRATOS_EXPECT_FALSE(divider_one_negative.mIsSplit);
    KRATOS_EXPECT_EQ(divider_one_negative.mDivisionsNumber, 1);
    KRATOS_EXPECT_EQ(divider_one_negative.mSplitEdgesNumber, 0);

    nodal_distances[2] = 0.3;
    auto divider_one_positive = DivideTriangle2D3<Node>(p_elem->GetGeometry(), nodal_distances);
    divider_one_positive.GenerateDivision();

    KRATOS_EXPECT_FALSE(divider_one_positive.mIsSplit);
    KRATOS_EXPECT_EQ(divider_one_positive.mDivisionsNumber, 1);
    KRATOS_EXPECT_EQ(divider_one_positive.mSplitEdgesNumber, 0);
}

}  // namespace Kratos::Testing.
