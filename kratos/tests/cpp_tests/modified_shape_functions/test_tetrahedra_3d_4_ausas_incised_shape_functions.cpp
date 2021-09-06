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
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions.h"

namespace Kratos
{
    namespace Testing
    {

        KRATOS_TEST_CASE_IN_SUITE(AusasIncisedShapeFunctionsTetrahedra3D4Horizontal, KratosCoreFastSuite)
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

            // Set the elemental distances vector
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 4> distances_vector;
            distances_vector[0] = -1.0;
            distances_vector[1] = -1.0;
            distances_vector[2] = -1.0;
            distances_vector[3] =  1.0;
            array_1d<double, 6> edge_dist_extra_vector;
            edge_dist_extra_vector[0] = -1.0;
            edge_dist_extra_vector[1] = -1.0;
            edge_dist_extra_vector[2] = -1.0;
            edge_dist_extra_vector[3] =  0.5;
            edge_dist_extra_vector[4] = -1.0;
            edge_dist_extra_vector[5] =  0.5;

            // Call the modified shape functions calculator
            Tetrahedra3D4AusasIncisedShapeFunctions tetrahedra_ausas_shape_functions(p_geometry, distances_vector, edge_dist_extra_vector);
            Matrix positive_side_sh_func, negative_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
            Vector positive_side_weights, negative_side_weights;

            tetrahedra_ausas_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
                positive_side_sh_func,
                positive_side_sh_func_gradients,
                positive_side_weights,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface modified shape functions calculator
            Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
            Vector positive_interface_side_weights, negative_interface_side_weights;

            tetrahedra_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_interface_side_sh_func,
                positive_interface_side_sh_func_gradients,
                positive_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

            tetrahedra_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            const std::vector<double> ref_pos_sh_func = {1.0/8.0, 0.0, 1.0/8.0, 3.0/4.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_side_sh_func, 0), ref_pos_sh_func, tolerance);

            Matrix ref_neg_sh_func(3, 4);
            const std::array<double, 12> ref_neg_sh_func_array = {1.0/8.0, 1.0/2.0, 1.0/4.0, 1.0/8.0,
                                                                  3.0/8.0, 1.0/4.0, 1.0/4.0, 1.0/8.0,
                                                                  1.0/8.0, 1.0/4.0, 3.0/8.0, 1.0/4.0};
            for (unsigned int i = 0; i < ref_neg_sh_func.size1(); i++) {
                for (unsigned int j = 0; j < ref_neg_sh_func.size2(); j++) {
                    ref_neg_sh_func(i, j) = ref_neg_sh_func_array[i * ref_neg_sh_func.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func, ref_neg_sh_func, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 0.02083333333, tolerance);
            const std::vector<double> ref_neg_weights = {0.04166666666, 0.08333333333, 0.02083333333};
            KRATOS_CHECK_VECTOR_NEAR(negative_side_weights, ref_neg_weights, tolerance);

            // Check Gauss pts. shape functions gradients values
            Matrix ref_pos_sh_func_grad_0(4, 3);
            Matrix ref_neg_sh_func_grad_0(4, 3);
            Matrix ref_neg_sh_func_grad_1(4, 3);
            Matrix ref_neg_sh_func_grad_2(4, 3);
            const std::array<double, 12> ref_pos_sh_func_grad_0_array = {-1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_sh_func_grad_0_array = {-1.0, -1.0, -1.0,  2.0,  1.0,  2.0,  0.0,  1.0,  0.0, -1.0, -1.0, -1.0};
            const std::array<double, 12> ref_neg_sh_func_grad_1_array = {-1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_sh_func_grad_2_array = {-1.0, -1.0, -1.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0, -1.0,  0.0,  1.0};
            for (unsigned int i = 0; i < ref_pos_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_sh_func_grad_0.size2(); j++) {
                    ref_pos_sh_func_grad_0(i, j) = ref_pos_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_0(i, j) = ref_neg_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_1(i, j) = ref_neg_sh_func_grad_1_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_2(i, j) = ref_neg_sh_func_grad_2_array[i * ref_pos_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func_gradients[0], ref_pos_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[0], ref_neg_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[1], ref_neg_sh_func_grad_1, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[2], ref_neg_sh_func_grad_2, tolerance);

            // Check interface shape function values
            const std::vector<double> ref_pos_interface_sh_func = {1.0/6.0,     0.0, 1.0/6.0, 2.0/3.0};
            const std::vector<double> ref_neg_interface_sh_func = {1.0/6.0, 1.0/3.0, 1.0/6.0, 1.0/3.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_interface_side_sh_func, 0), ref_pos_interface_sh_func, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(row(negative_interface_side_sh_func, 0), ref_neg_interface_sh_func, tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.125, tolerance);

            // Check Gauss pts. interface shape function gradients values
            Matrix ref_pos_interface_sh_func_grad_0(4, 3);
            Matrix ref_neg_interface_sh_func_grad_0(4, 3);
            const std::array<double, 12> ref_pos_interface_sh_func_grad_0_array = {-1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_interface_sh_func_grad_0_array = {-1.0, -1.0, -1.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0, -1.0,  0.0,  1.0};
            for (unsigned int i = 0; i < ref_pos_interface_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_interface_sh_func_grad_0.size2(); j++) {
                    ref_pos_interface_sh_func_grad_0(i, j) = ref_pos_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                    ref_neg_interface_sh_func_grad_0(i, j) = ref_neg_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_interface_side_sh_func_gradients[0], ref_pos_interface_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_interface_side_sh_func_gradients[0], ref_neg_interface_sh_func_grad_0, tolerance);

            // Check Gauss pts. outwards unit normal values
            const std::vector<double> ref_pos_area_normals = {0.0, 0.0, -0.125};
            const std::vector<double> ref_neg_area_normals = {0.0, 0.0,  0.125};
            KRATOS_CHECK_VECTOR_NEAR(positive_side_area_normals[0], ref_pos_area_normals, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(negative_side_area_normals[0], ref_neg_area_normals, tolerance);
        }


        KRATOS_TEST_CASE_IN_SUITE(AusasIncisedShapeFunctionsTetrahedra3D4Oblique, KratosCoreFastSuite)
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

            // Set the elemental distances vector
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 4> distances_vector;
            distances_vector[0] = -1.0;
            distances_vector[1] =  1.0;
            distances_vector[2] = -1.0;
            distances_vector[3] =  1.0;
            array_1d<double, 6> edge_dist_extra_vector;
            edge_dist_extra_vector[0] = -1.0;
            edge_dist_extra_vector[1] = -1.0;
            edge_dist_extra_vector[2] = -1.0;
            edge_dist_extra_vector[3] =  0.5;
            edge_dist_extra_vector[4] = -1.0;
            edge_dist_extra_vector[5] =  0.5;

            // Call the modified shape functions calculator
            Tetrahedra3D4AusasIncisedShapeFunctions tetrahedra_ausas_shape_functions(p_geometry, distances_vector, edge_dist_extra_vector);
            Matrix positive_side_sh_func, negative_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
            Vector positive_side_weights, negative_side_weights;

            tetrahedra_ausas_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
                positive_side_sh_func,
                positive_side_sh_func_gradients,
                positive_side_weights,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface modified shape functions calculator
            Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
            Vector positive_interface_side_weights, negative_interface_side_weights;

            tetrahedra_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_interface_side_sh_func,
                positive_interface_side_sh_func_gradients,
                positive_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<array_1d<double,3>> positive_side_area_normals, negative_side_area_normals;

            tetrahedra_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            Matrix ref_pos_sh_func(3, 4);
            Matrix ref_neg_sh_func(3, 4);
            const std::array<double, 12> ref_pos_sh_func_array = {1.0/8.0, 1.0/4.0, 1.0/8.0, 1.0/2.0,
                                                                      0.0, 3.0/4.0,     0.0, 1.0/4.0,
                                                                      0.0, 1.0/2.0, 1.0/8.0, 3.0/8.0};
            const std::array<double, 12> ref_neg_sh_func_array = {1.0/4.0,     0.0, 5.0/8.0, 1.0/8.0,
                                                                  3.0/8.0,     0.0, 3.0/8.0, 1.0/4.0,
                                                                  5.0/8.0,     0.0, 1.0/4.0, 1.0/8.0};
            for (unsigned int i = 0; i < ref_neg_sh_func.size1(); i++) {
                for (unsigned int j = 0; j < ref_neg_sh_func.size2(); j++) {
                    ref_pos_sh_func(i, j) = ref_pos_sh_func_array[i * ref_neg_sh_func.size2() + j];
                    ref_neg_sh_func(i, j) = ref_neg_sh_func_array[i * ref_neg_sh_func.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func, ref_pos_sh_func, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func, ref_neg_sh_func, tolerance);

            // Check Gauss pts. weights
            const std::vector<double> ref_pos_weights = {0.02083333333, 0.04166666666, 0.02083333333};
            const std::vector<double> ref_neg_weights = {0.02083333333, 0.02083333333, 0.04166666666};
            KRATOS_CHECK_VECTOR_NEAR(positive_side_weights, ref_pos_weights, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(negative_side_weights, ref_neg_weights, tolerance);

            // Check Gauss pts. shape functions gradients values
            Matrix ref_pos_sh_func_grad_0(4, 3);
            Matrix ref_pos_sh_func_grad_1(4, 3);
            Matrix ref_pos_sh_func_grad_2(4, 3);
            Matrix ref_neg_sh_func_grad_0(4, 3);
            Matrix ref_neg_sh_func_grad_1(4, 3);
            Matrix ref_neg_sh_func_grad_2(4, 3);
            const std::array<double, 12> ref_pos_sh_func_grad_0_array = {-2.0, -1.0, -1.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_pos_sh_func_grad_1_array = { 0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_pos_sh_func_grad_2_array = { 0.0,  0.0,  0.0,  2.0,  0.0,  0.0, -2.0,  0.0, -1.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_sh_func_grad_0_array = {-2.0, -2.0, -2.0,  0.0,  0.0,  0.0,  2.0,  2.0,  1.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_sh_func_grad_1_array = { 0.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_sh_func_grad_2_array = { 0.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0};
            for (unsigned int i = 0; i < ref_pos_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_sh_func_grad_0.size2(); j++) {
                    ref_pos_sh_func_grad_0(i, j) = ref_pos_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_pos_sh_func_grad_1(i, j) = ref_pos_sh_func_grad_1_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_pos_sh_func_grad_2(i, j) = ref_pos_sh_func_grad_2_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_0(i, j) = ref_neg_sh_func_grad_0_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_1(i, j) = ref_neg_sh_func_grad_1_array[i * ref_pos_sh_func_grad_0.size2() + j];
                    ref_neg_sh_func_grad_2(i, j) = ref_neg_sh_func_grad_2_array[i * ref_pos_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func_gradients[0], ref_pos_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func_gradients[1], ref_pos_sh_func_grad_1, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(positive_side_sh_func_gradients[2], ref_pos_sh_func_grad_2, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[0], ref_neg_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[1], ref_neg_sh_func_grad_1, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_side_sh_func_gradients[2], ref_neg_sh_func_grad_2, tolerance);

            // Check interface shape function values
            const std::vector<double> ref_pos_interface_sh_func_0 = {1.0/6.0, 1.0/3.0, 1.0/6.0, 1.0/3.0};
            const std::vector<double> ref_pos_interface_sh_func_1 = {    0.0, 2.0/3.0, 1.0/6.0, 1.0/6.0};
            const std::vector<double> ref_neg_interface_sh_func_0 = {1.0/3.0,     0.0, 1.0/2.0, 1.0/6.0};
            const std::vector<double> ref_neg_interface_sh_func_1 = {1.0/2.0,     0.0, 1.0/6.0, 1.0/3.0};
            KRATOS_CHECK_VECTOR_NEAR(row(positive_interface_side_sh_func, 0), ref_pos_interface_sh_func_0, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(row(positive_interface_side_sh_func, 1), ref_pos_interface_sh_func_1, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(row(negative_interface_side_sh_func, 0), ref_neg_interface_sh_func_0, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(row(negative_interface_side_sh_func, 1), ref_neg_interface_sh_func_1, tolerance);

            // Check interface Gauss pts. weights
            const std::vector<double> ref_interface_weights = {0.176777, 0.176777};
            KRATOS_CHECK_VECTOR_NEAR(positive_interface_side_weights, ref_interface_weights, 1e-6);
            KRATOS_CHECK_VECTOR_NEAR(negative_interface_side_weights, ref_interface_weights, 1e-6);

            // Check Gauss pts. interface shape function gradients values
            Matrix ref_pos_interface_sh_func_grad_0(4, 3);
            Matrix ref_pos_interface_sh_func_grad_1(4, 3);
            Matrix ref_neg_interface_sh_func_grad_0(4, 3);
            Matrix ref_neg_interface_sh_func_grad_1(4, 3);
            const std::array<double, 12> ref_pos_interface_sh_func_grad_0_array = {-2.0, -1.0, -1.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_pos_interface_sh_func_grad_1_array = { 0.0,  0.0,  0.0,  2.0,  0.0,  0.0, -2.0,  0.0, -1.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_interface_sh_func_grad_0_array = {-2.0, -2.0, -2.0,  0.0,  0.0,  0.0,  2.0,  2.0,  1.0,  0.0,  0.0,  1.0};
            const std::array<double, 12> ref_neg_interface_sh_func_grad_1_array = { 0.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0};
            for (unsigned int i = 0; i < ref_pos_interface_sh_func_grad_0.size1(); i++) {
                for (unsigned int j = 0; j < ref_pos_interface_sh_func_grad_0.size2(); j++) {
                    ref_pos_interface_sh_func_grad_0(i, j) = ref_pos_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                    ref_pos_interface_sh_func_grad_1(i, j) = ref_pos_interface_sh_func_grad_1_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                    ref_neg_interface_sh_func_grad_0(i, j) = ref_neg_interface_sh_func_grad_0_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                    ref_neg_interface_sh_func_grad_1(i, j) = ref_neg_interface_sh_func_grad_1_array[i * ref_pos_interface_sh_func_grad_0.size2() + j];
                }
            }
            KRATOS_CHECK_MATRIX_NEAR(positive_interface_side_sh_func_gradients[0], ref_pos_interface_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(positive_interface_side_sh_func_gradients[1], ref_pos_interface_sh_func_grad_1, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_interface_side_sh_func_gradients[0], ref_neg_interface_sh_func_grad_0, tolerance);
            KRATOS_CHECK_MATRIX_NEAR(negative_interface_side_sh_func_gradients[1], ref_neg_interface_sh_func_grad_1, tolerance);

            // Check Gauss pts. outwards unit normal values
            const std::vector<double> ref_pos_area_normals_0 = {-0.125, 0.0, -0.125};
            const std::vector<double> ref_pos_area_normals_1 = {-0.125, 0.0, -0.125};
            const std::vector<double> ref_neg_area_normals_0 = { 0.125, 0.0,  0.125};
            const std::vector<double> ref_neg_area_normals_1 = { 0.125, 0.0,  0.125};
            KRATOS_CHECK_VECTOR_NEAR(positive_side_area_normals[0], ref_pos_area_normals_0, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(positive_side_area_normals[1], ref_pos_area_normals_1, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(negative_side_area_normals[0], ref_neg_area_normals_0, tolerance);
            KRATOS_CHECK_VECTOR_NEAR(negative_side_area_normals[1], ref_neg_area_normals_1, tolerance);
        }


        KRATOS_TEST_CASE_IN_SUITE(AusasIncisedShapeFunctionsTetrahedra3D4Volumes, KratosCoreFastSuite)
        {
            Model current_model;
            // Generate a model part with the previous
            ModelPart& base_model_part = current_model.CreateModelPart("Tetrahedra");
            base_model_part.AddNodalSolutionStepVariable(DISTANCE);

            // Fill the model part geometry data
            base_model_part.CreateNewNode(1, 0.0, 2.0, 0.0);
            base_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            base_model_part.CreateNewNode(3, 1.0, 2.0, 0.0);
            base_model_part.CreateNewNode(4, 1.0, 2.0, 2.0);
            Properties::Pointer p_properties(new Properties(0));
            base_model_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

            // Set the elemental distances vector
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 4> distances_vector;
            distances_vector[0] = -0.5;
            distances_vector[1] =  1.0;
            distances_vector[2] = -0.5;
            distances_vector[3] =  1.0;
            array_1d<double, 6> edge_dist_extra_vector;
            edge_dist_extra_vector[0] = -1.0;
            edge_dist_extra_vector[1] = -1.0;
            edge_dist_extra_vector[2] = -1.0;
            edge_dist_extra_vector[3] =  1.0/3.0;
            edge_dist_extra_vector[4] = -1.0;
            edge_dist_extra_vector[5] =  1.0/3.0;

            // Call the modified shape functions calculator
            Tetrahedra3D4AusasIncisedShapeFunctions tetrahedra_ausas_shape_functions(p_geometry, distances_vector, edge_dist_extra_vector);
            Matrix positive_side_sh_func, negative_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, negative_side_sh_func_gradients;
            Vector positive_side_weights, negative_side_weights;

            tetrahedra_ausas_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
                positive_side_sh_func,
                positive_side_sh_func_gradients,
                positive_side_weights,
                GeometryData::GI_GAUSS_2);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::GI_GAUSS_2);

            const double tolerance = 1e-10;

            // Check Gauss pts. weights
            const unsigned int n_gauss_pos = positive_side_weights.size();
            const unsigned int n_gauss_neg = negative_side_weights.size();

            double pos_vol = 0.0;
            for (unsigned int i=0; i<n_gauss_pos; ++i) {
                pos_vol += positive_side_weights(i);
            }

            double neg_vol = 0.0;
            for (unsigned int i=0; i<n_gauss_neg; ++i) {
                neg_vol += negative_side_weights(i);
            }

            const double tot_vol = 2.0*1.0*2.0/6.0;
            KRATOS_CHECK_NEAR(pos_vol+neg_vol, tot_vol, tolerance);
        }
    }   // namespace Testing.
}  // namespace Kratos.
