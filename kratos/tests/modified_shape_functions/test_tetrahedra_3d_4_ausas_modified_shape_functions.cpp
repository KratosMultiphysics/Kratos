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
#include "includes/gid_io.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

namespace Kratos
{
    namespace Testing
    {

        KRATOS_TEST_CASE_IN_SUITE(AusasModifiedShapeFunctionsTetrahedra3D4Horizontal, KratosCoreFastSuite)
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
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 4> distances_vector;
            for (unsigned int i = 0; i < p_geometry->size(); ++i) {
                distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
            }

            base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

            Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

            // Call the modified shape functions calculator
            Tetrahedra3D4AusasModifiedShapeFunctions tetrahedra_ausas_shape_functions(p_geometry, r_elemental_distances);
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
            
            // Call the external face modified shape functions calculator
            Matrix pos_ext_face_sh_func_0, neg_ext_face_sh_func_0,
                   pos_ext_face_sh_func_1, neg_ext_face_sh_func_1,
                   pos_ext_face_sh_func_2, neg_ext_face_sh_func_2,
                   pos_ext_face_sh_func_3, neg_ext_face_sh_func_3;

            ModifiedShapeFunctions::ShapeFunctionsGradientsType 
                pos_ext_face_sh_func_gradients_0, neg_ext_face_sh_func_gradients_0,
                pos_ext_face_sh_func_gradients_1, neg_ext_face_sh_func_gradients_1,
                pos_ext_face_sh_func_gradients_2, neg_ext_face_sh_func_gradients_2,
                pos_ext_face_sh_func_gradients_3, neg_ext_face_sh_func_gradients_3;

            Vector pos_ext_face_weights_0, neg_ext_face_weights_0, 
                   pos_ext_face_weights_1, neg_ext_face_weights_1, 
                   pos_ext_face_weights_2, neg_ext_face_weights_2,
                   pos_ext_face_weights_3, neg_ext_face_weights_3;

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_0,
                pos_ext_face_sh_func_gradients_0,
                pos_ext_face_weights_0,
                0,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_0,
                neg_ext_face_sh_func_gradients_0,
                neg_ext_face_weights_0,
                0,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_1,
                pos_ext_face_sh_func_gradients_1,
                pos_ext_face_weights_1,
                1,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_1,
                neg_ext_face_sh_func_gradients_1,
                neg_ext_face_weights_1,
                1,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_2,
                pos_ext_face_sh_func_gradients_2,
                pos_ext_face_weights_2,
                2,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_2,
                neg_ext_face_sh_func_gradients_2,
                neg_ext_face_weights_2,
                2,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_3,
                pos_ext_face_sh_func_gradients_3,
                pos_ext_face_weights_3,
                3,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_3,
                neg_ext_face_sh_func_gradients_3,
                neg_ext_face_weights_3,
                3,
                GeometryData::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<Vector> positive_side_area_normals, negative_side_area_normals;

            tetrahedra_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::GI_GAUSS_1);

            // Call the exterior faces outwards normal area vector calculator
            std::vector<Vector>
                area_normals_pos_face_0, area_normals_neg_face_0,
                area_normals_pos_face_1, area_normals_neg_face_1,
                area_normals_pos_face_2, area_normals_neg_face_2,
                area_normals_pos_face_3, area_normals_neg_face_3;

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_0, 0, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_0, 0, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_1, 1, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_1, 1, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_2, 2, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_2, 2, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_3, 3, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_3, 3, GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,3), 1.00, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,3), 0.00, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,3), 0.00, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(2,0), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,1), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,2), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,3), 0.00, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 0.02083333333, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(0), 0.08333333333, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(1), 0.04166666666, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(2), 0.02083333333, tolerance);

            // Check Gauss pts. shape functions gradients values
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,2), 0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,2),  0.0, tolerance);

            // Check interface shape function values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,3), 1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,3), 0.0,     tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.125, tolerance);

            // Check Gauss pts. interface shape function gradients values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,2), 0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,2),  0.0, tolerance);

            // Check face 0 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,3), 1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](3,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](3,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](3,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_0(0), 0.216506, 10e-5);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,0),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,1), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(1,0),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(1,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(1,2), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(1,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[1](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_0(0), 0.433013, 10e-5);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_0(1), 0.216506, 10e-5);

            KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](2), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](0), 0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](1), 0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](2), 0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[1](0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[1](1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[1](2), 0.125, tolerance);

            // Check face 1 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,3), 1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](3,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](3,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](3,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_1(0), 0.125, tolerance);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,0), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,1),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,1),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,2), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_1(0),  0.25, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_1(1), 0.125, tolerance);

            KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](0), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](0), -0.25, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[1](0), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[1](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[1](2), 0.0, tolerance);

            // Check face 2 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,3), 1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](3,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](3,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](3,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_2(0), 0.125, tolerance);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,0), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,2),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(1,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(1,1), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(1,2),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(1,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](1,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](1,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[1](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_2(0),  0.25, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_2(1), 0.125, tolerance);

            KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](1), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](1), -0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[1](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[1](1), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[1](2), 0.0, tolerance);
            
            // Check face 3 values
            KRATOS_CHECK_EQUAL(pos_ext_face_weights_3.size(), 0);
            KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_3.size1(), 0);
            KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_3.size2(), 4);
            KRATOS_CHECK_EQUAL(pos_ext_face_sh_func_gradients_3.size(), 0);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_3(0), 0.5, tolerance);

            KRATOS_CHECK_EQUAL(area_normals_pos_face_3.size(), 0);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[0](2), -0.5, tolerance);

            // Check Gauss pts. outwards unit normal values
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), -0.125, tolerance);

            KRATOS_CHECK_NEAR(negative_side_area_normals[0](0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](2),  0.125, tolerance);
        }


        KRATOS_TEST_CASE_IN_SUITE(AusasModifiedShapeFunctionsTetrahedra3D4Oblique, KratosCoreFastSuite)
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
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 4> distances_vector;
            for (unsigned int i = 0; i < p_geometry->size(); ++i) {
                distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
            }

            base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

            Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

            // Call the modified shape functions calculator
            Tetrahedra3D4AusasModifiedShapeFunctions tetrahedra_ausas_shape_functions(p_geometry, r_elemental_distances);
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

            tetrahedra_ausas_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                negative_interface_side_sh_func,
                negative_interface_side_sh_func_gradients,
                negative_interface_side_weights,
                GeometryData::GI_GAUSS_1);

            // Call the external face modified shape functions calculator
            Matrix pos_ext_face_sh_func_0, neg_ext_face_sh_func_0,
                   pos_ext_face_sh_func_1, neg_ext_face_sh_func_1,
                   pos_ext_face_sh_func_2, neg_ext_face_sh_func_2,
                   pos_ext_face_sh_func_3, neg_ext_face_sh_func_3;

            ModifiedShapeFunctions::ShapeFunctionsGradientsType 
                pos_ext_face_sh_func_gradients_0, neg_ext_face_sh_func_gradients_0,
                pos_ext_face_sh_func_gradients_1, neg_ext_face_sh_func_gradients_1,
                pos_ext_face_sh_func_gradients_2, neg_ext_face_sh_func_gradients_2,
                pos_ext_face_sh_func_gradients_3, neg_ext_face_sh_func_gradients_3;

            Vector pos_ext_face_weights_0, neg_ext_face_weights_0, 
                   pos_ext_face_weights_1, neg_ext_face_weights_1, 
                   pos_ext_face_weights_2, neg_ext_face_weights_2,
                   pos_ext_face_weights_3, neg_ext_face_weights_3;

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_0,
                pos_ext_face_sh_func_gradients_0,
                pos_ext_face_weights_0,
                0,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_0,
                neg_ext_face_sh_func_gradients_0,
                neg_ext_face_weights_0,
                0,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_1,
                pos_ext_face_sh_func_gradients_1,
                pos_ext_face_weights_1,
                1,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_1,
                neg_ext_face_sh_func_gradients_1,
                neg_ext_face_weights_1,
                1,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_2,
                pos_ext_face_sh_func_gradients_2,
                pos_ext_face_weights_2,
                2,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_2,
                neg_ext_face_sh_func_gradients_2,
                neg_ext_face_weights_2,
                2,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                pos_ext_face_sh_func_3,
                pos_ext_face_sh_func_gradients_3,
                pos_ext_face_weights_3,
                3,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                neg_ext_face_sh_func_3,
                neg_ext_face_sh_func_gradients_3,
                neg_ext_face_weights_3,
                3,
                GeometryData::GI_GAUSS_1);

            // Call the interface outwards normal unit vector calculator
            std::vector<Vector> positive_side_area_normals, negative_side_area_normals;

            tetrahedra_ausas_shape_functions.ComputePositiveSideInterfaceAreaNormals(
                positive_side_area_normals,
                GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeSideInterfaceAreaNormals(
                negative_side_area_normals,
                GeometryData::GI_GAUSS_1);

            // Call the exterior faces outwards normal area vector calculator
            std::vector<Vector>
                area_normals_pos_face_0, area_normals_neg_face_0,
                area_normals_pos_face_1, area_normals_neg_face_1,
                area_normals_pos_face_2, area_normals_neg_face_2,
                area_normals_pos_face_3, area_normals_neg_face_3;

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_0, 0, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_0, 0, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_1, 1, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_1, 1, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_2, 2, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_2, 2, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputePositiveExteriorFaceAreaNormals(
                area_normals_pos_face_3, 3, GeometryData::GI_GAUSS_1);

            tetrahedra_ausas_shape_functions.ComputeNegativeExteriorFaceAreaNormals(
                area_normals_neg_face_3, 3, GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 0.25, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,3), 0.75, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,0), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,1), 0.50, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,2), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,3), 0.50, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,0), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,1), 0.75, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,2), 0.00, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,3), 0.25, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 0.75, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 0.00, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,3), 0.00, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 0.00, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 0.50, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,3), 0.00, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,0), 0.25, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,1), 0.00, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,2), 0.75, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,3), 0.00, tolerance);

            // Check Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_side_weights(0), 0.02083333333, tolerance);
            KRATOS_CHECK_NEAR(positive_side_weights(1), 0.02083333333, tolerance);
            KRATOS_CHECK_NEAR(positive_side_weights(2), 0.04166666666, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(0), 0.04166666666, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(1), 0.02083333333, tolerance);
            KRATOS_CHECK_NEAR(negative_side_weights(2), 0.02083333333, tolerance);

            // Check Gauss pts. shape functions gradients values
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](1,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,2),  0.0, tolerance);

            // Check interface shape function values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,3), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,0), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,1), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,2), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,3), 1.0/3.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,3), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,0), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,1), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,2), 0.0,     tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,3), 1.0/3.0, tolerance);

            // Check interface Gauss pts. weights
            KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.176777, 1e-6);
            KRATOS_CHECK_NEAR(positive_interface_side_weights(1), 0.176777, 1e-6);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.176777, 1e-6);
            KRATOS_CHECK_NEAR(negative_interface_side_weights(1), 0.176777, 1e-6);

            // Check Gauss pts. interface shape function gradients values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,2),  0.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](3,2),  0.0, tolerance);

            // Check face 0 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,0),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,2),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(0,3), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(1,0),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(1,1), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(1,2),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_0(1,3), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](1,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_0[1](3,2),  1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_0(0), 0.216506, 10e-5);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_0(1), 0.433013, 10e-5);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,2), 1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_0(0,3), 0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_0[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_0(0), 0.216506, 10e-5);

            KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[0](2), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[1](0), 0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[1](1), 0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_0[1](2), 0.250, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_0[0](2), 0.125, tolerance);

            // Check face 1 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_1(0,3), 1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_1[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_1(0), 0.125, tolerance);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,0), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,1),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(0,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,1),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,2), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_1(1,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](2,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_1[1](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_1(0),  0.25, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_1(1), 0.125, tolerance);

            KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](0), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_1[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](0), -0.25, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[1](0), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[1](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_1[1](2), 0.0, tolerance);

            // Check face 2 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,0),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,2),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(0,3), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(1,0),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(1,1), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(1,2),     0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_2(1,3), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](3,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](1,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_2[1](3,2),  1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_2(0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_2(1), 0.250, tolerance);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,0), 1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_2(0,3), 0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_2[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_2(0), 0.125, tolerance);

            KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](1), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[0](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[1](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[1](1), -0.25, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_2[1](2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](1), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_2[0](2), 0.0, tolerance);

            // Check face 3 values
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_3(0,0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_3(0,1), 1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_3(0,2), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_3(0,3), 0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](0,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](1,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_sh_func_gradients_3[0](3,2),  1.0, tolerance);
            KRATOS_CHECK_NEAR(pos_ext_face_weights_3(0), 0.125, tolerance);
            
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,0), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,1),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(0,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(1,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(1,1),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(1,2), 2.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_3(1,3),     0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](0,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[0](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](0,0), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](0,1), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](0,2), -2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](1,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](2,0),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](2,1),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](2,2),  2.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_sh_func_gradients_3[1](3,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_3(0),  0.25, tolerance);
            KRATOS_CHECK_NEAR(neg_ext_face_weights_3(1), 0.125, tolerance);

            KRATOS_CHECK_NEAR(area_normals_pos_face_3[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_3[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_pos_face_3[0](2), -0.125, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[0](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[0](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[0](2), -0.25, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[1](0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[1](1), 0.0, tolerance);
            KRATOS_CHECK_NEAR(area_normals_neg_face_3[1](2), -0.125, tolerance);

            // Check Gauss pts. outwards unit normal values
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](0), -0.125, 1e-6);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](1),  0.0, 1e-6);
            KRATOS_CHECK_NEAR(positive_side_area_normals[0](2), -0.125, 1e-6);
            KRATOS_CHECK_NEAR(positive_side_area_normals[1](0), -0.125, 1e-6);
            KRATOS_CHECK_NEAR(positive_side_area_normals[1](1),  0.0, 1e-6);
            KRATOS_CHECK_NEAR(positive_side_area_normals[1](2), -0.125, 1e-6);

            KRATOS_CHECK_NEAR(negative_side_area_normals[0](0), 0.125, 1e-6);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](1), 0.0, 1e-6);
            KRATOS_CHECK_NEAR(negative_side_area_normals[0](2), 0.125, 1e-6);
            KRATOS_CHECK_NEAR(negative_side_area_normals[1](0), 0.125, 1e-6);
            KRATOS_CHECK_NEAR(negative_side_area_normals[1](1), 0.0, 1e-6);
            KRATOS_CHECK_NEAR(negative_side_area_normals[1](2), 0.125, 1e-6);
        }


        KRATOS_TEST_CASE_IN_SUITE(AusasModifiedShapeFunctionsTetrahedra3D4Volumes, KratosCoreFastSuite)
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

            // Set the DISTANCE field
            base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = -0.5;
            base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) =  1.0;
            base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = -0.5;
            base_model_part.Nodes()[4].FastGetSolutionStepValue(DISTANCE) =  1.0;

            // Set the elemental distances vector
            Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();

            array_1d<double, 4> distances_vector;
            for (unsigned int i = 0; i < p_geometry->size(); ++i) {
                distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
            }

            base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

            Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

            // Call the modified shape functions calculator
            Tetrahedra3D4AusasModifiedShapeFunctions tetrahedra_ausas_shape_functions(p_geometry, r_elemental_distances);
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
