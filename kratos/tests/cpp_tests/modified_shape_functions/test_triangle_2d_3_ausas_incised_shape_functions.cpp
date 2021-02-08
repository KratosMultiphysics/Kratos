//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//  			 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
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
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

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
                GeometryData::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface modified shape functions calculator
            Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
            Vector positive_interface_side_weights, negative_interface_side_weights;

            triangle_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_interface_side_sh_func,
                positive_interface_side_sh_func_gradients,
                positive_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<Vector> positive_side_area_normals, negative_side_area_normals;

            triangle_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 1.0/6.0, tolerance);
              KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 5.0/6.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/6.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 0.0, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/8.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(1), 1.0/4.0, tolerance);

            // Check Gauss pts. shape functions gradients values
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  2.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  0.0, tolerance);

            // Check interface shape function values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 1.0/4.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 3.0/4.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 1.0/2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/4.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 1.0/4.0, tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.50, tolerance);

            // Check Gauss pts. interface shape function gradients values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  2.0, tolerance);

            // Check Gauss pts. outwards unit normal values
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](1), -0.5, tolerance);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), 0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_area_normals[0](0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](1),  0.5, tolerance);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](2),  0.0, tolerance);
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
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

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
                GeometryData::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
                negative_side_sh_func,
                negative_side_sh_func_gradients,
                negative_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface modified shape functions calculator
            Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
            Vector positive_interface_side_weights, negative_interface_side_weights;

            triangle_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_interface_side_sh_func,
                positive_interface_side_sh_func_gradients,
                positive_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<Vector> positive_side_area_normals, negative_side_area_normals;

            triangle_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::GI_GAUSS_1);

            triangle_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 5.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 1.0/6.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/2.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 1.0/3.0, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/8.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(1), 1.0/4.0, tolerance);

            // Check Gauss pts. shape functions gradients values
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, tolerance);

            // Check interface shape function values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 3.0/4.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 1.0/4.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 1.0/2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/4.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 1.0/4.0, tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.50, tolerance);

            // Check Gauss pts. interface shape function gradients values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);

            // Check Gauss pts. outwards unit normal values
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](0), -0.5, tolerance);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), 0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_area_normals[0](0),  0.5, tolerance);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](2),  0.0, tolerance);
        }
    }   // namespace Testing.
}  // namespace Kratos.
