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
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{
	namespace Testing
	{

		KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3Horizontal, KratosCoreFastSuite)
		{
			// Generate a model part with the previous
			ModelPart base_model_part("Triangle");
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
			Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();
			
			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < p_geometry->size(); ++i) {
				distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
			}
			
			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);
			
			const Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);
			
			// Call the modified shape functions calculator
			Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
			Matrix positive_side_sh_func, negative_side_sh_func;
			std::vector<Matrix> positive_side_sh_func_gradients, negative_side_sh_func_gradients;
			Vector positive_side_weights, negative_side_weights;

			triangle_shape_functions.GetPositiveSideShapeFunctionsAndGradientsValues(positive_side_sh_func,
																					 positive_side_sh_func_gradients,
																					 positive_side_weights,
																					 GeometryData::GI_GAUSS_1);

			triangle_shape_functions.GetNegativeSideShapeFunctionsAndGradientsValues(negative_side_sh_func,
																					 negative_side_sh_func_gradients,
																					 negative_side_weights,
																					 GeometryData::GI_GAUSS_1);
				
			// Check shape functions values
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 2.0/3.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/2.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/3.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 1.0/2.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 1.0/3.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 1.0/6.0, 1e-5);
			
			// Check Gauss pts. weights
			KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/8.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_weights(1), 1.0/4.0, 1e-5);
			
			// Check Gauss pts. shape functions gradients values
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, 1e-5);
				
		}
			
		KRATOS_TEST_CASE_IN_SUITE(ModifiedShapeFunctionsTriangle2D3Vertical, KratosCoreFastSuite)
		{
			// Generate a model part with the previous
			ModelPart base_model_part("Triangle");
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
			Geometry<Node<3>>::Pointer p_geometry = base_model_part.Elements()[1].pGetGeometry();
			
			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < p_geometry->size(); ++i) {
				distances_vector(i) = (*p_geometry)[i].FastGetSolutionStepValue(DISTANCE);
			}
			
			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);
			
			const Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);
			
			// Call the modified shape functions calculator
			Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
			Matrix positive_side_sh_func, negative_side_sh_func;
			std::vector<Matrix> positive_side_sh_func_gradients, negative_side_sh_func_gradients;
			Vector positive_side_weights, negative_side_weights;

			triangle_shape_functions.GetPositiveSideShapeFunctionsAndGradientsValues(positive_side_sh_func,
																					 positive_side_sh_func_gradients,
																					 positive_side_weights,
																					 GeometryData::GI_GAUSS_1);

			triangle_shape_functions.GetNegativeSideShapeFunctionsAndGradientsValues(negative_side_sh_func,																				
																				     negative_side_sh_func_gradients,
																					 negative_side_weights,
																					 GeometryData::GI_GAUSS_1);

			// Check shape functions values
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,0), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,1), 2.0/3.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func(0,2), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,0), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,1), 1.0/3.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(0,2), 1.0/2.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,0), 1.0/2.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,1), 1.0/6.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func(1,2), 1.0/3.0, 1e-5);

			// Check Gauss pts. weights
			KRATOS_CHECK_NEAR(positive_side_weights(0), 1.0/8.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_weights(0), 1.0/8.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_weights(1), 1.0/4.0, 1e-5);

			// Check Gauss pts. shape functions gradients values
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,0), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](0,1), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,0),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](1,1),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,0),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(positive_side_sh_func_gradients[0](2,1),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,0), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](0,1), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,0),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](1,1),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,0),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[0](2,1),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,0), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](0,1), -1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,0),  1.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](1,1),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,0),  0.0, 1e-5);
			KRATOS_CHECK_NEAR(negative_side_sh_func_gradients[1](2,1),  1.0, 1e-5);

			// Call the modified shape functions interface calculator

		
		}
        
	}   // namespace Testing.
}  // namespace Kratos.
