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
#include "utilities/discont_utils.h"

namespace Kratos
{
	namespace Testing
	{

		KRATOS_TEST_CASE_IN_SUITE(TriangleHorizontalDiscontUtils, KratosCoreFastSuite)
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

			// Compute the triangle intersection
			BoundedMatrix<double, 3, 2> point_coordinates;
			BoundedMatrix<double, 3, 2> continuous_N_gradients;
			array_1d<double, 3> nodal_distances;
			array_1d<double, 3> partition_volumes;
			BoundedMatrix<double, 3, 3> gauss_pt_continuous_N_values;
			array_1d<double, 3> partition_signs;
			std::vector<Matrix> enriched_N_gradients_values(3);
			BoundedMatrix<double, 3, 3> enriched_N_values;
			array_1d<double, 3> edge_areas;


			auto rGeom = base_model_part.Elements()[1].GetGeometry();
			for (unsigned int inode=0; inode<base_model_part.NumberOfNodes(); ++inode)
			{
				point_coordinates(inode, 0) = rGeom[inode].X();
				point_coordinates(inode, 1) = rGeom[inode].Y();

				nodal_distances(inode) = rGeom[inode].FastGetSolutionStepValue(DISTANCE);
			}

			unsigned int ndivisions = DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(point_coordinates,
																												 continuous_N_gradients,
																												 nodal_distances,
																												 partition_volumes,
																												 gauss_pt_continuous_N_values,
																												 partition_signs,
																												 enriched_N_gradients_values,
																												 enriched_N_values,
																												 edge_areas);

			// Number of divisions check
			KRATOS_CHECK_EQUAL(ndivisions, 3);

			// Continuous shape functions derivatives check
			KRATOS_CHECK_NEAR(continuous_N_gradients(0,0), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(0,1), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(1,0),  1.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(2,1),  1.0, 1e-6);

			// Partition volumes (areas) check
			KRATOS_CHECK_NEAR(partition_volumes(0), 0.250, 1e-6);
			KRATOS_CHECK_NEAR(partition_volumes(1), 0.125, 1e-6);
			KRATOS_CHECK_NEAR(partition_volumes(2), 0.125, 1e-6);
			const double total_volume = partition_volumes(0) + partition_volumes(1) + partition_volumes(2);
			KRATOS_CHECK_NEAR(total_volume, 0.5, 1e-6);

			// Gauss points shape function values check
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(0,0), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(0,1), 0.5, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(0,2), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(1,0), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(1,1), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(1,2), 2.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(2,0), 0.5, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(2,1), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(2,2), 1.0/3.0, 1e-6);

			// Check partition signs
			KRATOS_CHECK_EQUAL(partition_signs(0), -1);
			KRATOS_CHECK_EQUAL(partition_signs(1),  1);
			KRATOS_CHECK_EQUAL(partition_signs(2), -1);

			// Check partition gradients values
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](0,0), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](0,1), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](1,0),  1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](1,1),  1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](2,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](0,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](0,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](1,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](2,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](0,0), -2.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](0,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](1,0),  2.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](2,1),  0.0, 1e-6);

			// Check enriched shape function partition Gauss pts. values
			KRATOS_CHECK_NEAR(enriched_N_values(0,0), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(0,1), 2.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(0,2), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(1,0), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(1,1), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(1,2), 1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(2,0), 2.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(2,1), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(2,2), 0.0, 1e-6);

			// Check edge areas
			KRATOS_CHECK_NEAR(edge_areas(0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(edge_areas(1), 0.25, 1e-6);
			KRATOS_CHECK_NEAR(edge_areas(2), 0.25, 1e-6);
		}

		KRATOS_TEST_CASE_IN_SUITE(TriangleVerticalDiscontUtils, KratosCoreFastSuite)
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

			// Compute the triangle intersection
			BoundedMatrix<double, 3, 2> point_coordinates;
			BoundedMatrix<double, 3, 2> continuous_N_gradients;
			array_1d<double, 3> nodal_distances;
			array_1d<double, 3> partition_volumes;
			BoundedMatrix<double, 3, 3> gauss_pt_continuous_N_values;
			array_1d<double, 3> partition_signs;
			std::vector<Matrix> enriched_N_gradients_values(3);
			BoundedMatrix<double, 3, 3> enriched_N_values;
			array_1d<double, 3> edge_areas;


			auto rGeom = base_model_part.Elements()[1].GetGeometry();
			for (unsigned int inode=0; inode<base_model_part.NumberOfNodes(); ++inode)
			{
				point_coordinates(inode, 0) = rGeom[inode].X();
				point_coordinates(inode, 1) = rGeom[inode].Y();

				nodal_distances(inode) = rGeom[inode].FastGetSolutionStepValue(DISTANCE);
			}

			unsigned int ndivisions = DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(point_coordinates,
																												 continuous_N_gradients,
																												 nodal_distances,
																												 partition_volumes,
																												 gauss_pt_continuous_N_values,
																												 partition_signs,
																												 enriched_N_gradients_values,
																												 enriched_N_values,
																												 edge_areas);

			// Number of divisions check
			KRATOS_CHECK_EQUAL(ndivisions, 3);

			// Continuous shape functions derivatives check
			KRATOS_CHECK_NEAR(continuous_N_gradients(0,0), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(0,1), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(1,0),  1.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(continuous_N_gradients(2,1),  1.0, 1e-6);

			// Partition volumes (areas) check
			KRATOS_CHECK_NEAR(partition_volumes(0), 0.250, 1e-6);
			KRATOS_CHECK_NEAR(partition_volumes(1), 0.125, 1e-6);
			KRATOS_CHECK_NEAR(partition_volumes(2), 0.125, 1e-6);
			const double total_volume = partition_volumes(0) + partition_volumes(1) + partition_volumes(2);
			KRATOS_CHECK_NEAR(total_volume, 0.5, 1e-6);

			// Gauss points shape function values check
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(0,0), 0.5, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(0,1), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(0,2), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(1,0), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(1,1), 2.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(1,2), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(2,0), 1.0/6.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(2,1), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(gauss_pt_continuous_N_values(2,2), 0.5, 1e-6);

			// Check partition signs
			KRATOS_CHECK_EQUAL(partition_signs(0), -1);
			KRATOS_CHECK_EQUAL(partition_signs(1),  1);
			KRATOS_CHECK_EQUAL(partition_signs(2), -1);

			// Check partition gradients values
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](0,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](0,1), -1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](1,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[0](2,1),  1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](0,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](0,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](1,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](2,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[1](2,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](0,0), -2.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](0,1), -2.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](1,0),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](1,1),  0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](2,0),  2.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_gradients_values[2](2,1),  2.0, 1e-6);

			// Check enriched shape function partition Gauss pts. values
			KRATOS_CHECK_NEAR(enriched_N_values(0,0), 2.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(0,1), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(0,2), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(1,0), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(1,1), 1.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(1,2), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(2,0), 1.0/3.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(2,1), 0.0, 1e-6);
			KRATOS_CHECK_NEAR(enriched_N_values(2,2), 2.0/3.0, 1e-6);

			// Check edge areas
			KRATOS_CHECK_NEAR(edge_areas(0), 0.25, 1e-6);
			KRATOS_CHECK_NEAR(edge_areas(1), 0.25, 1e-6);
			KRATOS_CHECK_NEAR(edge_areas(2),  0.0, 1e-6);
		}

		KRATOS_TEST_CASE_IN_SUITE(TriangleNoIntersectionDiscontUtils, KratosCoreFastSuite)
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

			// Compute the triangle intersection
			BoundedMatrix<double, 3, 2> point_coordinates;
			BoundedMatrix<double, 3, 2> continuous_N_gradients;
			array_1d<double, 3> nodal_distances;
			array_1d<double, 3> partition_volumes;
			BoundedMatrix<double, 3, 3> gauss_pt_continuous_N_values;
			array_1d<double, 3> partition_signs;
			std::vector<Matrix> enriched_N_gradients_values(3);
			BoundedMatrix<double, 3, 3> enriched_N_values;
			array_1d<double, 3> edge_areas;


			auto rGeom = base_model_part.Elements()[1].GetGeometry();
			for (unsigned int inode=0; inode<base_model_part.NumberOfNodes(); ++inode)
			{
				point_coordinates(inode, 0) = rGeom[inode].X();
				point_coordinates(inode, 1) = rGeom[inode].Y();

				nodal_distances(inode) = rGeom[inode].FastGetSolutionStepValue(DISTANCE);
			}

			unsigned int ndivisions = DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(point_coordinates,
																												 continuous_N_gradients,
																												 nodal_distances,
																												 partition_volumes,
																												 gauss_pt_continuous_N_values,
																												 partition_signs,
																												 enriched_N_gradients_values,
																												 enriched_N_values,
																												 edge_areas);

			// Number of divisions check
			KRATOS_CHECK_EQUAL(ndivisions, 1);

			// Partition volumes (areas) check
			KRATOS_CHECK_NEAR(partition_volumes(0), 0.5, 1e-6);
		}
	}
}  // namespace Kratos.
