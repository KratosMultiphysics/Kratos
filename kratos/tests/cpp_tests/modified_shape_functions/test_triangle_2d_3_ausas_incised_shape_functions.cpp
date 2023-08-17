//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions.h"

namespace Kratos
{
    namespace Testing
    {

        KRATOS_TEST_CASE_IN_SUITE(AusasIncisedShapeFunctionsTriangle2D3Horizontal, KratosCoreFastSuite)
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

            // Set the elemental distances vector
            Geometry<Node>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 3> distances_vector;
            distances_vector[0] = -1.0;
            distances_vector[1] = -1.0;
            distances_vector[2] =  1.0;
            array_1d<double, 3> edge_dist_extra_vector;
            edge_dist_extra_vector[0] =  0.5;
            edge_dist_extra_vector[1] = -1.0;
            edge_dist_extra_vector[2] = -1.0;

            // Call the modified shape functions calculator
            Triangle2D3AusasIncisedShapeFunctions triangle_ausas_shape_functions(p_geometry, distances_vector, edge_dist_extra_vector);
            Matrix positive_side_sh_func, negative_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
            Vector positive_side_weights, negative_side_weights;

            triangle_ausas_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
                positive_side_sh_func,
                positive_side_sh_func_gradients,
                positive_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            // Call the interface modified shape functions calculator
            Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
            Vector positive_interface_side_weights, negative_interface_side_weights;

            triangle_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_interface_side_sh_func,
                positive_interface_side_sh_func_gradients,
                positive_interface_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

            triangle_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            const std::vector<double> ref_pos_sh_func = {0.0, 1.0/6.0, 5.0/6.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_side_sh_func, 0), ref_pos_sh_func, tolerance);

            Matrix ref_neg_sh_func(2, 3);
            const std::array<double, 6> ref_neg_sh_func_array = {1.0/3.0, 1.0/2.0, 1.0/6.0,
                                                                 2.0/3.0, 1.0/3.0,     0.0};
            for (unsigned int i = 0; i < ref_neg_sh_func.size1(); i++) {
                for (unsigned int j = 0; j < ref_neg_sh_func.size2(); j++) {
                    ref_neg_sh_func(i, j) = ref_neg_sh_func_array[i * ref_neg_sh_func.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func, ref_neg_sh_func, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, tolerance);
            const std::vector<double> ref_neg_weights = {1.0/8.0, 1.0/4.0};
            KRATOS_CHECK_VECTOR_NEAR(negative_side_weights, ref_neg_weights, tolerance);

            // Check Gauss pts. shape functions gradients values
            Matrix ref_pos_sh_func_grad_0(3, 2);
            Matrix ref_neg_sh_func_grad_0(3, 2);
            Matrix ref_neg_sh_func_grad_1(3, 2);
            const std::array<double, 6> ref_pos_sh_func_grad_0_array = { 0.0,  0.0,  1.0,  0.0, -1.0,  0.0};
            const std::array<double, 6> ref_neg_sh_func_grad_0_array = {-2.0, -2.0,  1.0,  0.0,  1.0,  2.0};
            const std::array<double, 6> ref_neg_sh_func_grad_1_array = {-1.0,  0.0,  1.0,  0.0,  0.0,  0.0};
            for (unsigned int i = 0; i < ref_pos_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_sh_func_grad_0.size2(); j++) {
                    ref_pos_sh_func_grad_0(i, j) = ref_pos_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_0(i, j) = ref_neg_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_1(i, j) = ref_neg_sh_func_grad_1_array[i * ref_pos_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func_gradients[0], ref_pos_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[0], ref_neg_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[1], ref_neg_sh_func_grad_1, tolerance);

            // Check interface shape function values
            const std::vector<double> ref_pos_interface_sh_func = {    0.0, 1.0/4.0, 3.0/4.0};
            const std::vector<double> ref_neg_interface_sh_func = {1.0/2.0, 1.0/4.0, 1.0/4.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_interface_side_sh_func, 0), ref_pos_interface_sh_func, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(row(negative_interface_side_sh_func, 0), ref_neg_interface_sh_func, tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.50, tolerance);

            // Check Gauss pts. interface shape function gradients values
            Matrix ref_pos_interface_sh_func_grad_0(3, 2);
            Matrix ref_neg_interface_sh_func_grad_0(3, 2);
            const std::array<double, 6> ref_pos_interface_sh_func_grad_0_array = { 0.0,  0.0,  1.0,  0.0, -1.0,  0.0};
            const std::array<double, 6> ref_neg_interface_sh_func_grad_0_array = {-2.0, -2.0,  1.0,  0.0,  1.0,  2.0};
            for (unsigned int i = 0; i < ref_pos_interface_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_interface_sh_func_grad_0.size2(); j++) {
                    ref_pos_interface_sh_func_grad_0(i, j) = ref_pos_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                    ref_neg_interface_sh_func_grad_0(i, j) = ref_neg_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_interface_side_sh_func_gradients[0], ref_pos_interface_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_interface_side_sh_func_gradients[0], ref_neg_interface_sh_func_grad_0, tolerance);

            // Check Gauss pts. outwards unit normal values
            const std::vector<double> ref_pos_area_normals = {0.0, -0.5, 0.0};
            const std::vector<double> ref_neg_area_normals = {0.0,  0.5, 0.0};
            KRATOS_CHECK_VECTOR_NEAR(positive_side_area_normals[0], ref_pos_area_normals, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(negative_side_area_normals[0], ref_neg_area_normals, tolerance);
        }


        KRATOS_TEST_CASE_IN_SUITE(AusasIncisedShapeFunctionsTriangle2D3Vertical, KratosCoreFastSuite)
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

            // Set the elemental distances vector
            Geometry<Node>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 3> distances_vector;
            distances_vector[0] = -1.0;
            distances_vector[1] =  1.0;
            distances_vector[2] = -1.0;
            array_1d<double, 3> edge_dist_extra_vector;
            edge_dist_extra_vector[0] =  0.5;
            edge_dist_extra_vector[1] = -1.0;
            edge_dist_extra_vector[2] = -1.0;

            // Call the modified shape functions calculator
            Triangle2D3AusasIncisedShapeFunctions triangle_ausas_shape_functions(p_geometry, distances_vector, edge_dist_extra_vector);
            Matrix positive_side_sh_func, negative_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
            Vector positive_side_weights, negative_side_weights;

            triangle_ausas_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
                positive_side_sh_func,
                positive_side_sh_func_gradients,
                positive_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            // Call the interface modified shape functions calculator
            Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
            Vector positive_interface_side_weights, negative_interface_side_weights;

            triangle_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_interface_side_sh_func,
                positive_interface_side_sh_func_gradients,
                positive_interface_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

            triangle_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::IntegrationMethod::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            const std::vector<double> ref_pos_sh_func = {0.0, 5.0/6.0, 1.0/6.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_side_sh_func, 0), ref_pos_sh_func, tolerance);

            Matrix ref_neg_sh_func(2, 3);
            const std::array<double, 6> ref_neg_sh_func_array = {1.0/3.0, 1.0/6.0, 1.0/2.0,
                                                                 2.0/3.0,     0.0, 1.0/3.0};
            for (unsigned int i = 0; i < ref_neg_sh_func.size1(); i++) {
                for (unsigned int j = 0; j < ref_neg_sh_func.size2(); j++) {
                    ref_neg_sh_func(i, j) = ref_neg_sh_func_array[i * ref_neg_sh_func.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func, ref_neg_sh_func, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, tolerance);
            const std::vector<double> ref_neg_weights = {1.0/8.0, 1.0/4.0};
            KRATOS_CHECK_VECTOR_NEAR(negative_side_weights, ref_neg_weights, tolerance);

            // Check Gauss pts. shape functions gradients values
            Matrix ref_pos_sh_func_grad_0(3, 2);
            Matrix ref_neg_sh_func_grad_0(3, 2);
            Matrix ref_neg_sh_func_grad_1(3, 2);
            const std::array<double, 6> ref_pos_sh_func_grad_0_array = { 0.0,  0.0,  0.0, -1.0,  0.0,  1.0};
            const std::array<double, 6> ref_neg_sh_func_grad_0_array = {-2.0, -2.0,  2.0,  1.0,  0.0,  1.0};
            const std::array<double, 6> ref_neg_sh_func_grad_1_array = { 0.0, -1.0,  0.0,  0.0,  0.0,  1.0};
            for (unsigned int i = 0; i < ref_pos_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_sh_func_grad_0.size2(); j++) {
                    ref_pos_sh_func_grad_0(i, j) = ref_pos_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_0(i, j) = ref_neg_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_1(i, j) = ref_neg_sh_func_grad_1_array[i * ref_pos_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func_gradients[0], ref_pos_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[0], ref_neg_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[1], ref_neg_sh_func_grad_1, tolerance);

            // Check interface shape function values
            const std::vector<double> ref_pos_interface_sh_func = {    0.0, 3.0/4.0, 1.0/4.0};
            const std::vector<double> ref_neg_interface_sh_func = {1.0/2.0, 1.0/4.0, 1.0/4.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_interface_side_sh_func, 0), ref_pos_interface_sh_func, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(row(negative_interface_side_sh_func, 0), ref_neg_interface_sh_func, tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.50, tolerance);

            // Check Gauss pts. interface shape function gradients values
            Matrix ref_pos_interface_sh_func_grad_0(3, 2);
            Matrix ref_neg_interface_sh_func_grad_0(3, 2);
            const std::array<double, 6> ref_pos_interface_sh_func_grad_0_array = { 0.0,  0.0,  0.0, -1.0,  0.0,  1.0};
            const std::array<double, 6> ref_neg_interface_sh_func_grad_0_array = {-2.0, -2.0,  2.0,  1.0,  0.0,  1.0};
            for (unsigned int i = 0; i < ref_pos_interface_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_interface_sh_func_grad_0.size2(); j++) {
                    ref_pos_interface_sh_func_grad_0(i, j) = ref_pos_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                    ref_neg_interface_sh_func_grad_0(i, j) = ref_neg_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_interface_side_sh_func_gradients[0], ref_pos_interface_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_interface_side_sh_func_gradients[0], ref_neg_interface_sh_func_grad_0, tolerance);

            // Check Gauss pts. outwards unit normal values
            const std::vector<double> ref_pos_area_normals = {-0.5, 0.0, 0.0};
            const std::vector<double> ref_neg_area_normals = { 0.5, 0.0, 0.0};
            KRATOS_CHECK_VECTOR_NEAR(positive_side_area_normals[0], ref_pos_area_normals, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(negative_side_area_normals[0], ref_neg_area_normals, tolerance);
        }
    }   // namespace Testing.
}  // namespace Kratos.
