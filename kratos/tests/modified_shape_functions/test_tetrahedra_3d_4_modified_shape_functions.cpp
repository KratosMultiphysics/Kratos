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
#include "includes/checks.h"
#include "includes/gid_io.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{
	namespace Testing
	{

		KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTetrahedra3D4Horizontal, KratosCoreFastSuite)
		{
			// Generate a model part with the previous
			ModelPart base_model_part("Tetrahedra");
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
			Tetrahedra3D4ModifiedShapeFunctions tetrahedra_shape_functions(p_geometry, r_elemental_distances);
			Matrix positive_side_sh_func, negative_side_sh_func;
			std::vector<Matrix> positive_side_sh_func_gradients, negative_side_sh_func_gradients;
			Vector positive_side_weights, negative_side_weights;

			tetrahedra_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
				positive_side_sh_func,
				positive_side_sh_func_gradients,
				positive_side_weights,
				GeometryData::GI_GAUSS_1);

			tetrahedra_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
				negative_side_sh_func,
				negative_side_sh_func_gradients,
				negative_side_weights,
				GeometryData::GI_GAUSS_1);

			// Call the interface modified shape functions calculator
			Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
			std::vector<Matrix> positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
			Vector positive_interface_side_weights, negative_interface_side_weights;

			tetrahedra_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
				positive_interface_side_sh_func,
				positive_interface_side_sh_func_gradients,
				positive_interface_side_weights,
				GeometryData::GI_GAUSS_1);

			tetrahedra_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
				negative_interface_side_sh_func,
				negative_interface_side_sh_func_gradients,
				negative_interface_side_weights,
				GeometryData::GI_GAUSS_1);
																							   
			// Call the interface outwards normal unit vector calculator
			std::vector<Vector> positive_side_unit_normals, negative_side_unit_normals;
			
			tetrahedra_shape_functions.ComputePositiveSideInterfaceUnitNormals(
				positive_side_unit_normals,
				GeometryData::GI_GAUSS_1);
			
			tetrahedra_shape_functions.ComputeNegativeSideInterfaceUnitNormals(
				negative_side_unit_normals,
				GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

			// Check shape functions values
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 0.125, tolerance);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 0.125, tolerance);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 0.125, tolerance);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,3), 0.625, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 0.375, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 0.250, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 0.250, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,3), 0.125, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 0.125, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 0.375, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 0.250, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,3), 0.250, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(2,0), 0.125, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(2,1), 0.125, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(2,2), 0.375, tolerance);
			KRATOS_CHECK_NEAR(negative_side_sh_func(2,3), 0.375, tolerance);

			// Check Gauss pts. weights
			KRATOS_CHECK_NEAR(positive_side_weights(0), 0.02083333333, tolerance);
			KRATOS_CHECK_NEAR(negative_side_weights(0), 0.08333333333, tolerance);
			KRATOS_CHECK_NEAR(negative_side_weights(1), 0.04166666666, tolerance);
			KRATOS_CHECK_NEAR(negative_side_weights(2), 0.02083333333, tolerance);

            // Check Gauss pts. shape functions gradients values
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,2),  1.0, tolerance);
            
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,2),  1.0, tolerance);

			// Check interface shape function values
			KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 1.0/6.0, tolerance);
			KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 1.0/6.0, tolerance);
			KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 1.0/6.0, tolerance);
			KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,3),     0.5, tolerance);
			KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 1.0/6.0, tolerance);
			KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/6.0, tolerance);
			KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 1.0/6.0, tolerance);
			KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,3),     0.5, tolerance);

			// Check interface Gauss pts. weights
			KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.125, tolerance);
			KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.125, tolerance);

            // Check Gauss pts. interface shape function gradients values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,2),  1.0, tolerance);

			// Check Gauss pts. outwards unit normal values
			KRATOS_CHECK_NEAR(positive_side_unit_normals[0](0),  0.0, tolerance);
			KRATOS_CHECK_NEAR(positive_side_unit_normals[0](1),  0.0, tolerance);
			KRATOS_CHECK_NEAR(positive_side_unit_normals[0](2), -1.0, tolerance);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[0](0),  0.0, tolerance);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[0](1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_unit_normals[0](2),  1.0, tolerance);
		}


		KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTetrahedra3D4Oblique, KratosCoreFastSuite)
		{
			// Generate a model part with the previous
			ModelPart base_model_part("Tetrahedra");
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
			Tetrahedra3D4ModifiedShapeFunctions tetrahedra_shape_functions(p_geometry, r_elemental_distances);
			Matrix positive_side_sh_func, negative_side_sh_func;
			std::vector<Matrix> positive_side_sh_func_gradients, negative_side_sh_func_gradients;
			Vector positive_side_weights, negative_side_weights;

			tetrahedra_shape_functions.ComputePositiveSideShapeFunctionsAndGradientsValues(
				positive_side_sh_func,
				positive_side_sh_func_gradients,
				positive_side_weights,
				GeometryData::GI_GAUSS_1);

			tetrahedra_shape_functions.ComputeNegativeSideShapeFunctionsAndGradientsValues(
				negative_side_sh_func,
				negative_side_sh_func_gradients,
				negative_side_weights,
				GeometryData::GI_GAUSS_1);

			// Call the interface modified shape functions calculator
			Matrix positive_interface_side_sh_func, negative_interface_side_sh_func;
			std::vector<Matrix> positive_interface_side_sh_func_gradients, negative_interface_side_sh_func_gradients;
			Vector positive_interface_side_weights, negative_interface_side_weights;

			tetrahedra_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
				positive_interface_side_sh_func,
				positive_interface_side_sh_func_gradients,
				positive_interface_side_weights,
				GeometryData::GI_GAUSS_1);

			tetrahedra_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
				negative_interface_side_sh_func,
				negative_interface_side_sh_func_gradients,
				negative_interface_side_weights,
				GeometryData::GI_GAUSS_1);
																							   
			// Call the interface outwards normal unit vector calculator
			std::vector<Vector> positive_side_unit_normals, negative_side_unit_normals;
			
			tetrahedra_shape_functions.ComputePositiveSideInterfaceUnitNormals(
				positive_side_unit_normals,
				GeometryData::GI_GAUSS_1);
			
			tetrahedra_shape_functions.ComputeNegativeSideInterfaceUnitNormals(
				negative_side_unit_normals,
				GeometryData::GI_GAUSS_1);

            const double tolerance = 1e-10;

            // Check shape functions values
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 0.250, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 0.125, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(0,3), 0.500, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,1), 0.250, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,2), 0.250, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(1,3), 0.375, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,1), 0.500, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,2), 0.125, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func(2,3), 0.250, tolerance);
            
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 0.500, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 0.250, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(0,3), 0.125, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 0.250, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 0.125, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 0.375, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(1,3), 0.250, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,0), 0.125, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,1), 0.250, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,2), 0.500, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func(2,3), 0.125, tolerance);

			// Check Gauss pts. weights
			KRATOS_CHECK_NEAR(positive_side_weights(0), 0.02083333333, tolerance);
			KRATOS_CHECK_NEAR(positive_side_weights(1), 0.02083333333, tolerance);
			KRATOS_CHECK_NEAR(positive_side_weights(2), 0.04166666666, tolerance);
			KRATOS_CHECK_NEAR(negative_side_weights(0), 0.04166666666, tolerance);
			KRATOS_CHECK_NEAR(negative_side_weights(1), 0.02083333333, tolerance);
			KRATOS_CHECK_NEAR(negative_side_weights(2), 0.02083333333, tolerance);

			// Check Gauss pts. shape functions gradients values
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[1](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[2](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[2](3,2),  1.0, tolerance);

            // Check interface shape function values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,1), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,2), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(0,3), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,0), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func(1,3), 1.0/6.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,0), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,1), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,2), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(0,3), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,0), 1.0/6.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,1), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,2), 1.0/3.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func(1,3), 1.0/6.0, tolerance);

            // Check interface Gauss pts. weights
			KRATOS_CHECK_NEAR(positive_interface_side_weights(0), 0.176777, 1e-6);
			KRATOS_CHECK_NEAR(positive_interface_side_weights(1), 0.176777, 1e-6);
			KRATOS_CHECK_NEAR(negative_interface_side_weights(0), 0.176777, 1e-6);
			KRATOS_CHECK_NEAR(negative_interface_side_weights(1), 0.176777, 1e-6);

			// Check Gauss pts. interface shape function gradients values
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[0](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(positive_interface_side_sh_func_gradients[1](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[0](3,2),  1.0, tolerance);

            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](0,0), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](0,1), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](0,2), -1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](1,0),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](1,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](2,1),  1.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](2,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](3,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](3,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(negative_interface_side_sh_func_gradients[1](3,2),  1.0, tolerance);

			// Check Gauss pts. outwards unit normal values
			KRATOS_CHECK_NEAR(positive_side_unit_normals[0](0), -0.707107, 1e-6);
			KRATOS_CHECK_NEAR(positive_side_unit_normals[0](1),       0.0, 1e-6);
			KRATOS_CHECK_NEAR(positive_side_unit_normals[0](2), -0.707107, 1e-6);
			KRATOS_CHECK_NEAR(positive_side_unit_normals[1](0), -0.707107, 1e-6);
			KRATOS_CHECK_NEAR(positive_side_unit_normals[1](1),       0.0, 1e-6);
            KRATOS_CHECK_NEAR(positive_side_unit_normals[1](2), -0.707107, 1e-6);
            
			KRATOS_CHECK_NEAR(negative_side_unit_normals[0](0),  0.707107, 1e-6);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[0](1),       0.0, 1e-6);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[0](2),  0.707107, 1e-6);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[1](0),  0.707107, 1e-6);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[1](1),       0.0, 1e-6);
			KRATOS_CHECK_NEAR(negative_side_unit_normals[1](2),  0.707107, 1e-6);
		}

	}   // namespace Testing.
}  // namespace Kratos.
