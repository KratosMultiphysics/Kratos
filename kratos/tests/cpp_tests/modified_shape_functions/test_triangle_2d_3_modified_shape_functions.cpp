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
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos::Testing
{

Triangle2D3ModifiedShapeFunctions::UniquePointer SetTriangle2D3ModifiedShapeFunctionsPointer()
{
    // Generate a test model part
    Model current_model;
    ModelPart& base_model_part = current_model.CreateModelPart("Triangle");
    base_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Fill the model part geometry data
    base_model_part.CreateNewNode(1, 0.0, 2.0, 0.0);
    base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    base_model_part.CreateNewNode(3, 1.0, 2.0, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    base_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

    // Set the DISTANCE field
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = -0.5;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) =  1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = -0.5;

    // Set the elemental distances vector
    auto p_geometry = base_model_part.Elements()[1].pGetGeometry();
    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < p_geometry->size(); ++i) {
        distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
    }
    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);
    const Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Create and return the modified shape functions instance
    auto p_modified_shape_functions_triangle = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(p_geometry, r_elemental_distances);
    return p_modified_shape_functions_triangle;
}

KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3Horizontal, KratosCoreFastSuite)
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
    Geometry<Node>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < p_geometry->size(); ++i) {
        distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    const Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Call the modified shape functions calculator
    Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
    Matrix positive_side_sh_func, negative_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
    Vector positive_side_weights, negative_side_weights;

    triangle_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
        negative_side_sh_func,
        negative_side_sh_func_gradients,
        negative_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the interface modified shape functions calculator
    Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
    Vector positive_interface_side_weights, negative_interface_side_weights;

    triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        positive_interface_side_sh_func,
        positive_interface_side_sh_func_gradients,
        positive_interface_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        negative_interface_side_sh_func,
        negative_interface_side_sh_func_gradients,
        negative_interface_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the external face modified shape functions calculator
    Matrix pos_ext_face_sh_func_0, neg_ext_face_sh_func_0,
            pos_ext_face_sh_func_1, neg_ext_face_sh_func_1,
            pos_ext_face_sh_func_2, neg_ext_face_sh_func_2;

    ModifiedShapeFunctions::ShapeFunctionsGradientsType
        pos_ext_face_sh_func_gradients_0, neg_ext_face_sh_func_gradients_0,
        pos_ext_face_sh_func_gradients_1, neg_ext_face_sh_func_gradients_1,
        pos_ext_face_sh_func_gradients_2, neg_ext_face_sh_func_gradients_2;

    Vector pos_ext_face_weights_0, neg_ext_face_weights_0,
            pos_ext_face_weights_1, neg_ext_face_weights_1,
            pos_ext_face_weights_2, neg_ext_face_weights_2;

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_0,
        pos_ext_face_sh_func_gradients_0,
        pos_ext_face_weights_0,
        0,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_0,
        neg_ext_face_sh_func_gradients_0,
        neg_ext_face_weights_0,
        0,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_1,
        pos_ext_face_sh_func_gradients_1,
        pos_ext_face_weights_1,
        1,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_1,
        neg_ext_face_sh_func_gradients_1,
        neg_ext_face_weights_1,
        1,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_2,
        pos_ext_face_sh_func_gradients_2,
        pos_ext_face_weights_2,
        2,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_2,
        neg_ext_face_sh_func_gradients_2,
        neg_ext_face_weights_2,
        2,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the interface outwards normal area vector calculator
    std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

    triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(
        positive_side_area_normals,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
        negative_side_area_normals,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the exterior faces outwards normal area vector calculator
    std::vector<array_1d<double,3>>
        area_normals_pos_face_0, area_normals_neg_face_0,
        area_normals_pos_face_1, area_normals_neg_face_1,
        area_normals_pos_face_2, area_normals_neg_face_2;

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_0, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_0, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_1, 1, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_1, 1, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_2, 2, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_2, 2, GeometryData::IntegrationMethod::GI_GAUSS_1);

    const double tolerance = 1e-10;

    // Check shape functions values
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 2.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 1.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 1.0/6.0, tolerance);

    // Check Gauss pts. weights
    KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/8.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_weights(1), 1.0/4.0, tolerance);

    // Check Gauss pts. shape functions gradients values
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, tolerance);

    // Check interface shape function values
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 0.25, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 0.25, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 0.50, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 0.25, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 0.25, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 0.50, tolerance);

    // Check interface Gauss pts. weights
    KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.50, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.50, tolerance);

    // Check Gauss pts. interface shape function gradients values
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);

    // Check Gauss pts. outwards area normal values
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](0), 0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](1), -0.5, tolerance);
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](0), 0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](2), 0.0, tolerance);

    // Check face 0 values
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,1), 0.25, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,2), 0.75, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_weights_0(0), 0.5*std::sqrt(2.0), tolerance);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,1), 0.75, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,2), 0.25, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_0(0), 0.5*std::sqrt(2.0), tolerance);

    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](2), 0.0, tolerance);

    // Check face 1 values
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,0), 0.25, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,2), 0.75, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_weights_1(0), 0.5, tolerance);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,0), 0.75, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,2), 0.25, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_1(0), 0.5, tolerance);

    KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](0), -0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](2),  0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](0), -0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](2),  0.0, tolerance);

    // Check face 2 values
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_2.size1(), 0);
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_2.size2(), 3);
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_gradients_2.size(), 0);
    KRATOS_CHECK_EQUAL(pos_ext_face_weights_2.size(), 0);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_2(0), 1.0, tolerance);

    KRATOS_CHECK_EQUAL(area_normals_pos_face_2.size(), 0);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](2),  0.0, tolerance);
}


KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3Vertical, KratosCoreFastSuite)
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
    Geometry<Node>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < p_geometry->size(); ++i) {
        distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    const Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Call the modified shape functions calculator
    Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
    Matrix positive_side_sh_func, negative_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
    Vector positive_side_weights, negative_side_weights;

    triangle_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
        negative_side_sh_func,
        negative_side_sh_func_gradients,
        negative_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the interface modified shape functions calculator
    Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
    Vector positive_interface_side_weights, negative_interface_side_weights;

    triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        positive_interface_side_sh_func,
        positive_interface_side_sh_func_gradients,
        positive_interface_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        negative_interface_side_sh_func,
        negative_interface_side_sh_func_gradients,
        negative_interface_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the interface outwards normal unit vector calculator
    std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

    triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(
        positive_side_area_normals,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
        negative_side_area_normals,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the external face modified shape functions calculator
    Matrix pos_ext_face_sh_func_0, neg_ext_face_sh_func_0,
            pos_ext_face_sh_func_1, neg_ext_face_sh_func_1,
            pos_ext_face_sh_func_2, neg_ext_face_sh_func_2;

    ModifiedShapeFunctions::ShapeFunctionsGradientsType
        pos_ext_face_sh_func_gradients_0, neg_ext_face_sh_func_gradients_0,
        pos_ext_face_sh_func_gradients_1, neg_ext_face_sh_func_gradients_1,
        pos_ext_face_sh_func_gradients_2, neg_ext_face_sh_func_gradients_2;

    Vector pos_ext_face_weights_0, neg_ext_face_weights_0,
            pos_ext_face_weights_1, neg_ext_face_weights_1,
            pos_ext_face_weights_2, neg_ext_face_weights_2;

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_0,
        pos_ext_face_sh_func_gradients_0,
        pos_ext_face_weights_0,
        0,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_0,
        neg_ext_face_sh_func_gradients_0,
        neg_ext_face_weights_0,
        0,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_1,
        pos_ext_face_sh_func_gradients_1,
        pos_ext_face_weights_1,
        1,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_1,
        neg_ext_face_sh_func_gradients_1,
        neg_ext_face_weights_1,
        1,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_2,
        pos_ext_face_sh_func_gradients_2,
        pos_ext_face_weights_2,
        2,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_2,
        neg_ext_face_sh_func_gradients_2,
        neg_ext_face_weights_2,
        2,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the exterior faces outwards normal area vector calculator
    std::vector<array_1d<double,3>>
        area_normals_pos_face_0, area_normals_neg_face_0,
        area_normals_pos_face_1, area_normals_neg_face_1,
        area_normals_pos_face_2, area_normals_neg_face_2;

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_0, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_0, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_1, 1, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_1, 1, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_2, 2, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_2, 2, GeometryData::IntegrationMethod::GI_GAUSS_1);

    const double tolerance = 1e-10;

    // Check shape functions values
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 2.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 1.0/3.0, tolerance);

    // Check Gauss pts. weights
    KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/8.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_weights(1), 1.0/4.0, tolerance);

    // Check Gauss pts. shape functions gradients values
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, tolerance);

    // Check interface shape function values
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 0.25, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 0.50, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 0.25, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 0.25, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 0.50, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 0.25, tolerance);

    // Check interface Gauss pts. weights
    KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.50, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.50, tolerance);

    // Check Gauss pts. interface shape function gradients values
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);

    // Check Gauss pts. outwards area normal values
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](0), -0.5, tolerance);
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](1), 0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](1), 0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](2), 0.0, tolerance);

    // Check face 0 values
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,1), 0.75, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,2), 0.25, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_weights_0(0), 0.5*std::sqrt(2.0), tolerance);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,1), 0.25, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,2), 0.75, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_0(0), 0.5*std::sqrt(2.0), tolerance);

    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](2), 0.0, tolerance);

    // Check face 1 values
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_1.size1(), 0);
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_1.size2(), 3);
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_gradients_1.size(), 0);
    KRATOS_CHECK_EQUAL(pos_ext_face_weights_1.size(), 0);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,1), 0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,2), 0.5, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_1(0), 1.0, tolerance);

    KRATOS_CHECK_EQUAL(area_normals_pos_face_1.size(), 0);
    KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](1), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](2), 0.0, tolerance);

    // Check face 2 values
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,0), 0.25, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,1), 0.75, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,2),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_weights_2(0), 0.5, tolerance);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,0), 0.75, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,1), 0.25, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,2),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_2(0), 0.5, tolerance);

    KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](0), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](1), -0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](0), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](1), -0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](2), 0.0, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3ZeroNode, KratosCoreFastSuite)
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
    base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) =  0.0;
    base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
    base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

    // Set the elemental distances vector
    Geometry<Node>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

    array_1d<double, 3> distances_vector;
    for (unsigned int i = 0; i < p_geometry->size(); ++i) {
        distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
    }

    base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

    const Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

    // Call the modified shape functions calculator
    Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
    Matrix positive_side_sh_func, negative_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
    Vector positive_side_weights, negative_side_weights;

    triangle_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
        negative_side_sh_func,
        negative_side_sh_func_gradients,
        negative_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the interface modified shape functions calculator
    Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
    Vector positive_interface_side_weights, negative_interface_side_weights;

    triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        positive_interface_side_sh_func,
        positive_interface_side_sh_func_gradients,
        positive_interface_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        negative_interface_side_sh_func,
        negative_interface_side_sh_func_gradients,
        negative_interface_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the external face modified shape functions calculator
    Matrix pos_ext_face_sh_func_0, neg_ext_face_sh_func_0,
            pos_ext_face_sh_func_1, neg_ext_face_sh_func_1,
            pos_ext_face_sh_func_2, neg_ext_face_sh_func_2;

    ModifiedShapeFunctions::ShapeFunctionsGradientsType
        pos_ext_face_sh_func_gradients_0, neg_ext_face_sh_func_gradients_0,
        pos_ext_face_sh_func_gradients_1, neg_ext_face_sh_func_gradients_1,
        pos_ext_face_sh_func_gradients_2, neg_ext_face_sh_func_gradients_2;

    Vector pos_ext_face_weights_0, neg_ext_face_weights_0,
            pos_ext_face_weights_1, neg_ext_face_weights_1,
            pos_ext_face_weights_2, neg_ext_face_weights_2;

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_0,
        pos_ext_face_sh_func_gradients_0,
        pos_ext_face_weights_0,
        0,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_0,
        neg_ext_face_sh_func_gradients_0,
        neg_ext_face_weights_0,
        0,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_1,
        pos_ext_face_sh_func_gradients_1,
        pos_ext_face_weights_1,
        1,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_1,
        neg_ext_face_sh_func_gradients_1,
        neg_ext_face_weights_1,
        1,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        pos_ext_face_sh_func_2,
        pos_ext_face_sh_func_gradients_2,
        pos_ext_face_weights_2,
        2,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        neg_ext_face_sh_func_2,
        neg_ext_face_sh_func_gradients_2,
        neg_ext_face_weights_2,
        2,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the interface outwards normal area vector calculator
    std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

    triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(
        positive_side_area_normals,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
        negative_side_area_normals,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    // Call the exterior faces outwards normal area vector calculator
    std::vector<array_1d<double,3>>
        area_normals_pos_face_0, area_normals_neg_face_0,
        area_normals_pos_face_1, area_normals_neg_face_1,
        area_normals_pos_face_2, area_normals_neg_face_2;

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_0, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_0, 0, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_1, 1, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_1, 1, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputePositiveExteriorFaceAreaNormals(
        area_normals_pos_face_2, 2, GeometryData::IntegrationMethod::GI_GAUSS_1);

    triangle_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
        area_normals_neg_face_2, 2, GeometryData::IntegrationMethod::GI_GAUSS_1);

    const double tolerance = 1e-10;

    KRATOS_CHECK_EQUAL(positive_side_sh_func.size1(), 1);
    KRATOS_CHECK_EQUAL(negative_side_sh_func.size1(), 1);

    // Check shape functions values
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 1.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 1.0/6.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/3.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/6.0, tolerance);

    // Check Gauss pts. weights
    KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/4.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/4.0, tolerance);

    // Check Gauss pts. shape functions gradients values
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);

    // Check interface shape function values
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 1.0/4.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 1.0/4.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 1.0/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/4.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 1.0/4.0, tolerance);

    // Check interface Gauss pts. weights
    KRATOS_CHECK_NEAR(positive_interface_side_weights(0), std::sqrt(2.0)/2.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_weights(0), std::sqrt(2.0)/2.0, tolerance);

    // Check Gauss pts. interface shape function gradients values
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);

    // Check Gauss pts. outwards area normal values
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](1), -0.5, tolerance);
    KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](0), -0.5, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(negative_side_area_normals[0](2), 0.0, tolerance);

    // Check face 0 values
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,1), 0.25, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,2), 0.75, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_weights_0(0), 0.5*std::sqrt(2.0), tolerance);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,1), 0.75, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,2), 0.25, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_0(0), 0.5*std::sqrt(2.0), tolerance);

    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](2), 0.0, tolerance);

    // Check face 1 values
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,1), 0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,2), 0.5, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(pos_ext_face_weights_1(0), 1.0, tolerance);

    KRATOS_CHECK_EQUAL(neg_ext_face_sh_func_1.size1(), 0);
    KRATOS_CHECK_EQUAL(neg_ext_face_sh_func_1.size2(), 3);
    KRATOS_CHECK_EQUAL(neg_ext_face_sh_func_gradients_1.size(), 0);
    KRATOS_CHECK_EQUAL(neg_ext_face_weights_1.size(), 0);

    KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](2),  0.0, tolerance);
    KRATOS_CHECK_EQUAL(area_normals_neg_face_1.size(), 0);

    // Check face 2 values
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_2.size1(), 0);
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_2.size2(), 3);
    KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_gradients_2.size(), 0);
    KRATOS_CHECK_EQUAL(pos_ext_face_weights_2.size(), 0);

    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,0), 0.5, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,1), 0.5, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,2), 0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,0),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,1),  1.0, tolerance);
    KRATOS_CHECK_NEAR(neg_ext_face_weights_2(0), 1.0, tolerance);

    KRATOS_CHECK_EQUAL(area_normals_pos_face_2.size(), 0);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](0),  0.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](2),  0.0, tolerance);
}


KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3Areas, KratosCoreFastSuite)
{
    auto p_triangle_shape_functions = SetTriangle2D3ModifiedShapeFunctionsPointer();
    Matrix positive_side_sh_func, negative_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
    Vector positive_side_weights, negative_side_weights;

    p_triangle_shape_functions->ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    p_triangle_shape_functions->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        negative_side_sh_func,
        negative_side_sh_func_gradients,
        negative_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    const double tolerance = 1e-10;

    // Check Gauss pts. weights
    const unsigned int n_gauss_pos = positive_side_weights.size();
    const unsigned int n_gauss_neg = negative_side_weights.size();

    double pos_area = 0.0;
    for (unsigned int i=0; i<n_gauss_pos; ++i) {
        pos_area += positive_side_weights(i);
    }

    double neg_area = 0.0;
    for (unsigned int i=0; i<n_gauss_neg; ++i) {
        neg_area += negative_side_weights(i);
    }

    const double tot_area = 2.0*1.0/2.0;
    KRATOS_CHECK_NEAR(pos_area+neg_area, tot_area, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3ComputeShapeFunctionsAndWeights, KratosCoreFastSuite)
{
    auto p_triangle_shape_functions = SetTriangle2D3ModifiedShapeFunctionsPointer();
    Matrix positive_side_sh_func, negative_side_sh_func;
    Vector positive_side_weights, negative_side_weights;

    p_triangle_shape_functions->ComputePositiveSideShapeFunctionsAndWeights(
        positive_side_sh_func,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    p_triangle_shape_functions->ComputeNegativeSideShapeFunctionsAndWeights(
        negative_side_sh_func,
        negative_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    const double tolerance = 1e-10;
    Matrix ref_pos_N(3,3);
    ref_pos_N(0,0) = 0.111111111111; ref_pos_N(0,1) = 0.4444444444444; ref_pos_N(0,2) = 0.444444444444;
    ref_pos_N(1,0) = 0.444444444444; ref_pos_N(1,1) = 0.4444444444444; ref_pos_N(1,2) = 0.111111111111;
    ref_pos_N(2,0) = 0.111111111111; ref_pos_N(2,1) = 0.7777777777777; ref_pos_N(2,2) = 0.111111111111;
    Matrix ref_neg_N(6,3);
    ref_neg_N(0,0) = 0.111111111111; ref_neg_N(0,1) = 0.2777777777777; ref_neg_N(0,2) = 0.611111111111;
    ref_neg_N(1,0) = 0.111111111111; ref_neg_N(1,1) = 0.1111111111111; ref_neg_N(1,2) = 0.777777777777;
    ref_neg_N(2,0) = 0.444444444444; ref_neg_N(2,1) = 0.2777777777777; ref_neg_N(2,2) = 0.277777777777;
    ref_neg_N(3,0) = 0.611111111111; ref_neg_N(3,1) = 0.2222222222222; ref_neg_N(3,2) = 0.166666666666;
    ref_neg_N(4,0) = 0.277777777777; ref_neg_N(4,1) = 0.0555555555555; ref_neg_N(4,2) = 0.666666666666;
    ref_neg_N(5,0) = 0.777777777777; ref_neg_N(5,1) = 0.0555555555555; ref_neg_N(5,2) = 0.166666666666;

    Vector ref_pos_weigt(3);
    ref_pos_weigt(0) = 0.148148148148;
    ref_pos_weigt(1) = 0.148148148148;
    ref_pos_weigt(2) = 0.148148148148;
    Vector ref_neg_weigt(6);
    ref_neg_weigt(0) = 0.074074074074;
    ref_neg_weigt(1) = 0.074074074074;
    ref_neg_weigt(2) = 0.074074074074;
    ref_neg_weigt(3) = 0.111111111111;
    ref_neg_weigt(4) = 0.111111111111;
    ref_neg_weigt(5) = 0.111111111111;
    KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func, ref_pos_N, tolerance);
    KRATOS_CHECK_VECTOR_NEAR(positive_side_weights, ref_pos_weigt, tolerance);
    KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func, ref_neg_N, tolerance);
    KRATOS_CHECK_VECTOR_NEAR(negative_side_weights, ref_neg_weigt, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3DomainSize, KratosCoreFastSuite)
{
    auto p_triangle_shape_functions = SetTriangle2D3ModifiedShapeFunctionsPointer();
    const double pos_dom_size = p_triangle_shape_functions->ComputePositiveSideDomainSize();
    const double neg_dom_size = p_triangle_shape_functions->ComputeNegativeSideDomainSize();

    const double tolerance = 1e-10;
    KRATOS_CHECK_NEAR(pos_dom_size, 0.444444444444, tolerance);
    KRATOS_CHECK_NEAR(neg_dom_size, 0.555555555555, tolerance);
    KRATOS_CHECK_NEAR(pos_dom_size+neg_dom_size, 1.0, tolerance);
}

}  // namespace Kratos::Testing.
