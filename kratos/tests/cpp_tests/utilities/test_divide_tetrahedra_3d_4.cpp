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
#include "includes/checks.h"
#include "utilities/divide_tetrahedra_3d_4.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4Horizontal, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Tetrahedra");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    base_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[4].FastGetSolutionStepValue(DISTANCE) =  1.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 4> distances_vector;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the tetrahedra splitting utility
    DivideTetrahedra3D4<Node> tetrahedra_splitter(r_geometry, r_elemental_distances);

    // Call the divide geometry method
    tetrahedra_splitter.GenerateDivision();

    // Call the intersection generation method
    tetrahedra_splitter.GenerateIntersectionsSkin();

    // Call the positive exterior faces generation method
    std::vector < unsigned int > pos_ext_faces_parent_ids;
    std::vector < DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType > pos_ext_faces;
    tetrahedra_splitter.GenerateExteriorFaces(
        pos_ext_faces,
        pos_ext_faces_parent_ids,
        tetrahedra_splitter.GetPositiveSubdivisions());

    // Call the negative exterior faces generation method
    std::vector < unsigned int > neg_ext_faces_parent_ids;
    std::vector < DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType > neg_ext_faces;
    tetrahedra_splitter.GenerateExteriorFaces(
        neg_ext_faces,
        neg_ext_faces_parent_ids,
        tetrahedra_splitter.GetNegativeSubdivisions());

    const double tolerance = 1e-10;

    // Check general splitting values
    KRATOS_CHECK(tetrahedra_splitter.mIsSplit);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 4);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 3);

    // Check split edges
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[0],  0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[1],  1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[2],  2);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[3],  3);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[4], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[5], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[6],  6);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[7], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[8],  8);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[9],  9);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[10],-1);

    // Check subdivisions
    const auto& r_positive_subdivision_0 = *(tetrahedra_splitter.GetPositiveSubdivisions()[0]);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].Z(), 0.5, tolerance);

    const auto& r_negative_subdivision_0 = *(tetrahedra_splitter.GetNegativeSubdivisions()[0]);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].Z(), 0.5, tolerance);

    // Check interfaces
    const auto& r_positive_interface_0 = *(tetrahedra_splitter.GetPositiveInterfaces()[0]);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].Z(), 0.5, tolerance);

    const auto& r_negative_interface_0 = *(tetrahedra_splitter.GetNegativeInterfaces()[0]);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].Z(), 0.5, tolerance);

    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetPositiveInterfacesParentIds()[0], 0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetNegativeInterfacesParentIds()[0], 2);

    // Check exterior faces
    KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 3);
    KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 7);

    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 0);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 0);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[2], 0);

    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 0);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 2);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 1);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[3], 2);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[4], 0);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[5], 1);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[6], 1);

    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Z(), 0.5, tolerance);

    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Z(), 0.0, tolerance);

}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4Oblique, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Tetrahedra");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    base_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) =  1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[4].FastGetSolutionStepValue(DISTANCE) =  1.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 4> distances_vector;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the tetrahedra splitting utility
    DivideTetrahedra3D4<Node> tetrahedra_splitter(r_geometry, r_elemental_distances);

    // Call the divide geometry method
    tetrahedra_splitter.GenerateDivision();

    // Call the intersection generation method
    tetrahedra_splitter.GenerateIntersectionsSkin();

    // Call the positive exterior faces generation method
    std::vector < unsigned int > pos_ext_faces_parent_ids;
    std::vector < DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType > pos_ext_faces;
    tetrahedra_splitter.GenerateExteriorFaces(
        pos_ext_faces,
        pos_ext_faces_parent_ids,
        tetrahedra_splitter.GetPositiveSubdivisions());

    // Call the negative exterior faces generation method
    std::vector < unsigned int > neg_ext_faces_parent_ids;
    std::vector < DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType > neg_ext_faces;
    tetrahedra_splitter.GenerateExteriorFaces(
        neg_ext_faces,
        neg_ext_faces_parent_ids,
        tetrahedra_splitter.GetNegativeSubdivisions());

    const double tolerance = 1e-10;

    // Check general splitting values
    KRATOS_CHECK(tetrahedra_splitter.mIsSplit);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 6);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 4);

    // Check split edges
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[0],  0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[1],  1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[2],  2);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[3],  3);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[4],  4);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[5], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[6],  6);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[7],  7);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[8], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[9],  9);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[10],-1);

    // Check subdivisions
    const auto &r_positive_subdivision_0 = *(tetrahedra_splitter.GetPositiveSubdivisions()[0]);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].Z(), 0.5, tolerance);

    const auto &r_negative_subdivision_0 = *(tetrahedra_splitter.GetNegativeSubdivisions()[0]);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].Z(), 0.0, tolerance);

    // Check interfaces
    const auto& r_positive_interface_0 = *(tetrahedra_splitter.GetPositiveInterfaces()[0]);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].Z(), 0.5, tolerance);

    const auto& r_positive_interface_1 = *(tetrahedra_splitter.GetPositiveInterfaces()[1]);
    KRATOS_CHECK_NEAR(r_positive_interface_1[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[2].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_1[2].Z(), 0.0, tolerance);

    const auto& r_negative_interface_0 = *(tetrahedra_splitter.GetNegativeInterfaces()[0]);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].Z(), 0.5, tolerance);

    const auto& r_negative_interface_1 = *(tetrahedra_splitter.GetNegativeInterfaces()[1]);
    KRATOS_CHECK_NEAR(r_negative_interface_1[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[1].Z(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_1[2].Z(), 0.5, tolerance);

    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetPositiveInterfacesParentIds()[0], 0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetPositiveInterfacesParentIds()[1], 2);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetNegativeInterfacesParentIds()[0], 0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetNegativeInterfacesParentIds()[1], 1);

    // Check exterior faces
    KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 6);
    KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 6);

    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 1);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 2);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[2], 0);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[3], 0);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[4], 1);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[5], 1);

    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 0);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 1);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 2);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[3], 2);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[4], 0);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[5], 2);

    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Z(), 0.5, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4NoDivision, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Tetrahedra");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    base_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;
    base_model_part.Nodes()[4].FastGetSolutionStepValue(DISTANCE) = 1.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 4> distances_vector;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the tetrahedra splitting utility
    DivideTetrahedra3D4<Node> tetrahedra_splitter(r_geometry, r_elemental_distances);

    // Call the divide geometry method
    tetrahedra_splitter.GenerateDivision();

    // Check general splitting values
    KRATOS_CHECK_IS_FALSE(tetrahedra_splitter.mIsSplit);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 0);

}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4ZeroNodes, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Tetrahedra");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    base_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) =  0.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;
    base_model_part.Nodes()[4].FastGetSolutionStepValue(DISTANCE) =  0.0;

    // Set the elemental distances vector
    Geometry < Node >& r_geometry = base_model_part.Elements()[1].GetGeometry();

    array_1d<double, 4> distances_vector;
    for (unsigned int i = 0; i < r_geometry.size(); ++i) {
        distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Build the tetrahedra splitting utility
    DivideTetrahedra3D4<Node> tetrahedra_splitter(r_geometry, r_elemental_distances);

    // Call the divide geometry method
    tetrahedra_splitter.GenerateDivision();

    // Call the intersection generation method
    tetrahedra_splitter.GenerateIntersectionsSkin();

    // Call the positive exterior faces generation method
    std::vector < unsigned int > pos_ext_faces_parent_ids;
    std::vector < DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType > pos_ext_faces;
    tetrahedra_splitter.GenerateExteriorFaces(
        pos_ext_faces,
        pos_ext_faces_parent_ids,
        tetrahedra_splitter.GetPositiveSubdivisions());

    // Call the negative exterior faces generation method
    std::vector < unsigned int > neg_ext_faces_parent_ids;
    std::vector < DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType > neg_ext_faces;
    tetrahedra_splitter.GenerateExteriorFaces(
        neg_ext_faces,
        neg_ext_faces_parent_ids,
        tetrahedra_splitter.GetNegativeSubdivisions());

    const double tolerance = 1e-10;

    // Check general splitting values
    KRATOS_CHECK(tetrahedra_splitter.mIsSplit);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 2);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 1);

    // Check split edges
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[0],  0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[1],  1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[2],  2);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[3],  3);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[4], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[5], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[6], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[7],  7);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[8], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[9], -1);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[10],-1);

    // Check subdivisions
    const auto& r_positive_subdivision_0 = *(tetrahedra_splitter.GetPositiveSubdivisions()[0]);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[2].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_subdivision_0[3].Z(), 0.0, tolerance);

    const auto& r_negative_subdivision_0 = *(tetrahedra_splitter.GetNegativeSubdivisions()[0]);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[2].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_subdivision_0[3].Z(), 0.0, tolerance);

    // Check interfaces
    const auto& r_positive_interface_0 = *(tetrahedra_splitter.GetPositiveInterfaces()[0]);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_positive_interface_0[2].Z(), 1.0, tolerance);

    const auto& r_negative_interface_0 = *(tetrahedra_splitter.GetNegativeInterfaces()[0]);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR(r_negative_interface_0[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetPositiveInterfacesParentIds()[0], 0);
    KRATOS_CHECK_EQUAL(tetrahedra_splitter.GetNegativeInterfacesParentIds()[0], 0);

    // Check exterior faces
    KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 3);
    KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 3);

    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 0);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 0);
    KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[2], 0);

    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 0);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 0);
    KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 0);

    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].Y(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Z(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].Z(), 0.0, tolerance);

    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].X(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Y(), 0.5, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].X(), 1.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Z(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].X(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].Y(), 0.0, tolerance);
    KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].Z(), 0.0, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4ConformantFaces, KratosCoreFastSuite)
{
    Model current_model;

    // Generate a model part with the previous
    ModelPart& base_model_part = current_model.CreateModelPart("Tetrahedra");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(195, -0.349832, 0.0522541, 0.808642);
    base_model_part.CreateNewNode(181, -0.132464, 0.144186, 0.788109);
    base_model_part.CreateNewNode(152, -0.232261, 0.0807726, 1.00046);
    base_model_part.CreateNewNode(170, -0.157697, -0.0702038, 0.799622);
    base_model_part.CreateNewNode(163, -0.335408, -0.1450987, 0.920880);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element3D4N", 471, {195, 181, 152, 170}, p_properties);
    base_model_part.CreateNewElement("Element3D4N", 526, {170, 195, 163, 152}, p_properties);

    // Set the DISTANCE field
    base_model_part.GetNode(195).FastGetSolutionStepValue(DISTANCE) = 0.181358;
    base_model_part.GetNode(181).FastGetSolutionStepValue(DISTANCE) = 0.201891;
    base_model_part.GetNode(152).FastGetSolutionStepValue(DISTANCE) = -0.0104574;
    base_model_part.GetNode(170).FastGetSolutionStepValue(DISTANCE) = 0.190378;
    base_model_part.GetNode(163).FastGetSolutionStepValue(DISTANCE) = 0.06912;

    // Set the elemental distances vector
    for (auto &i_elem : base_model_part.Elements()) {
        auto &r_geometry = i_elem.GetGeometry();
        array_1d<double, 4> distances_vector;
        for (unsigned int i = 0; i < r_geometry.size(); ++i) {
            distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
        }
        i_elem.SetValue(ELEMENTAL_DISTANCES, distances_vector);
    }

    // Build the tetrahedra splitting utilities
    auto p_elem_471 = base_model_part.pGetElement(471);
    auto p_elem_526 = base_model_part.pGetElement(526);

    DivideTetrahedra3D4<Node> tetra_split_471(p_elem_471->GetGeometry(), p_elem_471->GetValue(ELEMENTAL_DISTANCES));
    DivideTetrahedra3D4<Node> tetra_split_526(p_elem_526->GetGeometry(), p_elem_526->GetValue(ELEMENTAL_DISTANCES));

    // Call the divide geometry method
    tetra_split_471.GenerateDivision();
    tetra_split_526.GenerateDivision();

    // Call the positive exterior faces generation method
    std::vector<unsigned int> pos_ext_faces_parent_ids_471;
    std::vector<unsigned int> pos_ext_faces_parent_ids_526;
    std::vector<DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType> pos_ext_faces_471;
    std::vector<DivideTetrahedra3D4<Node>::IndexedPointGeometryPointerType> pos_ext_faces_526;
    tetra_split_471.GenerateExteriorFaces(pos_ext_faces_471, pos_ext_faces_parent_ids_471, tetra_split_471.GetPositiveSubdivisions());
    tetra_split_526.GenerateExteriorFaces(pos_ext_faces_526, pos_ext_faces_parent_ids_526, tetra_split_526.GetPositiveSubdivisions());

    // Check that shared positive faces are have the same splitting pattern
    KRATOS_CHECK_NEAR((*(pos_ext_faces_471[3])).Area(), (*(pos_ext_faces_526[4])).Area(), 1.0e-12);
    KRATOS_CHECK_NEAR((*(pos_ext_faces_471[2])).Area(), (*(pos_ext_faces_526[5])).Area(), 1.0e-12);
}

}  // namespace Kratos::Testing.
