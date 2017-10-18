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
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/geometry_splitting_utils.h"

namespace Kratos
{
	namespace Testing
	{

		KRATOS_TEST_CASE_IN_SUITE(TriangleHorizontalGeometrySplittingUtils, KratosCoreFastSuite)
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
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the triangle splitting utility
			TriangleSplittingUtils triangle_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			TriangleSplittingUtils::IndexedPointsContainerType aux_points_set;
			std::vector < TriangleSplittingUtils::IndexedPointGeometryType > positive_subdivisions;
			std::vector < TriangleSplittingUtils::IndexedPointGeometryType > negative_subdivisions;

			bool is_divided = triangle_splitter.DivideGeometry(aux_points_set, positive_subdivisions, negative_subdivisions);

			// Call the shape functions computation method
			Matrix positive_side_sh_func, negative_side_sh_func;
			Vector positive_side_weights, negative_side_weights; 
			triangle_splitter.GetShapeFunctionValues(positive_side_sh_func, positive_side_weights, positive_subdivisions, GeometryData::GI_GAUSS_1);
			triangle_splitter.GetShapeFunctionValues(negative_side_sh_func, negative_side_weights, negative_subdivisions, GeometryData::GI_GAUSS_1);

			// Check general splitting values
			KRATOS_CHECK(is_divided);
			KRATOS_CHECK_EQUAL(triangle_splitter.mDivisionsNumber, 3);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdgesNumber, 2);

			// Check split edges
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[0],  0);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[1],  1);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[2],  2);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[3], -1);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[4],  4);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[5],  5);

			// Check subdivisions
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[0].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[0].Y(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[1].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[1].Y(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[2].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[2].Y(), 1.0, 1e-5);

			KRATOS_CHECK_NEAR((negative_subdivisions[0])[0].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[0].Y(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[1].X(), 1.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[1].Y(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[2].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[2].Y(), 0.5, 1e-5);

			KRATOS_CHECK_NEAR((negative_subdivisions[1])[0].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[0].Y(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[1].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[1].Y(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[2].X(), 1.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[2].Y(), 0.0, 1e-5);

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
			
		}
		
		KRATOS_TEST_CASE_IN_SUITE(TriangleVerticalGeometrySplittingUtils, KratosCoreFastSuite)
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
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();
			
			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}
			
			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);
			
			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);
			
			// Build the triangle splitting utility
			TriangleSplittingUtils triangle_splitter(r_geometry, r_elemental_distances);
			
			// Call the divide geometry method
			TriangleSplittingUtils::IndexedPointsContainerType aux_points_set;
			std::vector < TriangleSplittingUtils::IndexedPointGeometryType > positive_subdivisions;
			std::vector < TriangleSplittingUtils::IndexedPointGeometryType > negative_subdivisions;
			
			bool is_divided = triangle_splitter.DivideGeometry(aux_points_set, positive_subdivisions, negative_subdivisions);
			
			// Call the shape functions computation method
			Matrix positive_side_sh_func, negative_side_sh_func;
			Vector positive_side_weights, negative_side_weights; 
			triangle_splitter.GetShapeFunctionValues(positive_side_sh_func, positive_side_weights, positive_subdivisions, GeometryData::GI_GAUSS_1);
			triangle_splitter.GetShapeFunctionValues(negative_side_sh_func, negative_side_weights, negative_subdivisions, GeometryData::GI_GAUSS_1);
			
			// Check general splitting values
			KRATOS_CHECK(is_divided);
			KRATOS_CHECK_EQUAL(triangle_splitter.mDivisionsNumber, 3);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdgesNumber, 2);
			
			// Check split edges
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[0],  0);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[1],  1);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[2],  2);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[3],  3);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[4],  4);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdges[5], -1);
			
			// Check subdivisions
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[0].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[0].Y(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[1].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[1].Y(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[2].X(), 1.0, 1e-5);
			KRATOS_CHECK_NEAR((positive_subdivisions[0])[2].Y(), 0.0, 1e-5);
			
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[0].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[0].Y(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[1].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[1].Y(), 1.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[2].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[0])[2].Y(), 0.0, 1e-5);
			
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[0].X(), 0.5, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[0].Y(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[1].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[1].Y(), 1.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[2].X(), 0.0, 1e-5);
			KRATOS_CHECK_NEAR((negative_subdivisions[1])[2].Y(), 0.0, 1e-5);
			
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

		}
		
		KRATOS_TEST_CASE_IN_SUITE(TriangleNoIntersectionGeometrySplittingUtils, KratosCoreFastSuite)
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
			base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
			base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
			base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

			// Set the elemental distances vector
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.PointsNumber(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the triangle splitting utility
			TriangleSplittingUtils triangle_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			TriangleSplittingUtils::IndexedPointsContainerType aux_points_set;
			std::vector < TriangleSplittingUtils::IndexedPointGeometryType > positive_subdivisions;
			std::vector < TriangleSplittingUtils::IndexedPointGeometryType > negative_subdivisions;

			bool is_divided = triangle_splitter.DivideGeometry(aux_points_set, positive_subdivisions, negative_subdivisions);

			// Check general splitting values
			KRATOS_CHECK_IS_FALSE(is_divided);
			KRATOS_CHECK_EQUAL(triangle_splitter.mDivisionsNumber, 1);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdgesNumber, 0);

		}
	}
}  // namespace Kratos.
