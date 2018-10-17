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

namespace Kratos
{
	namespace Testing
	{

		KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3Horizontal, KratosCoreFastSuite)
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
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the triangle splitting utility
			DivideTriangle2D3 triangle_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			triangle_splitter.GenerateDivision();

			// Call the intersection generation method
			triangle_splitter.GenerateIntersectionsSkin();

			// Call the positive exterior faces generation method
			std::vector < unsigned int > pos_ext_faces_parent_ids;
			std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > pos_ext_faces;
			triangle_splitter.GenerateExteriorFaces(
				pos_ext_faces,
				pos_ext_faces_parent_ids,
				triangle_splitter.mPositiveSubdivisions);

			// Call the negative exterior faces generation method
			std::vector < unsigned int > neg_ext_faces_parent_ids;
			std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > neg_ext_faces;
			triangle_splitter.GenerateExteriorFaces(
				neg_ext_faces,
				neg_ext_faces_parent_ids,
				triangle_splitter.mNegativeSubdivisions);

			const double tolerance = 1e-10;

			// Check general splitting values
			KRATOS_CHECK(triangle_splitter.mIsSplit);
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
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[2].Y(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[2].Y(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[2].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[2].Y(), 0.0, tolerance);

			// Check interfaces
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[1].Y(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[1].Y(), 0.5, tolerance);

			KRATOS_CHECK_EQUAL(triangle_splitter.mPositiveInterfacesParentIds[0], 0);
			KRATOS_CHECK_EQUAL(triangle_splitter.mNegativeInterfacesParentIds[0], 0);

			// Check exterior faces
			KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 2);
			KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 3);

			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 0);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 0);

			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 1);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 1);

			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Y(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Y(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Y(), 0.0, tolerance);
		}

		KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3Vertical, KratosCoreFastSuite)
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
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the triangle splitting utility
			DivideTriangle2D3 triangle_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			triangle_splitter.GenerateDivision();

			// Call the intersection generation method
			triangle_splitter.GenerateIntersectionsSkin();

			// Call the positive exterior faces generation method
			std::vector < unsigned int > pos_ext_faces_parent_ids;
			std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > pos_ext_faces;
			triangle_splitter.GenerateExteriorFaces(
				pos_ext_faces,
				pos_ext_faces_parent_ids,
				triangle_splitter.mPositiveSubdivisions);
		
			// Call the negative exterior faces generation method
			std::vector < unsigned int > neg_ext_faces_parent_ids;
			std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > neg_ext_faces;
			triangle_splitter.GenerateExteriorFaces(
				neg_ext_faces,
				neg_ext_faces_parent_ids,
				triangle_splitter.mNegativeSubdivisions);
		
			const double tolerance = 1e-10;

			// Check general splitting values
			KRATOS_CHECK(triangle_splitter.mIsSplit);
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
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[2].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveSubdivisions[0])[2].Y(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[0])[2].Y(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeSubdivisions[1])[2].Y(), 0.0, tolerance);

			// Check interfaces
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mPositiveInterfaces[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*triangle_splitter.mNegativeInterfaces[0])[1].Y(), 0.5, tolerance);

			KRATOS_CHECK_EQUAL(triangle_splitter.mPositiveInterfacesParentIds[0], 0);
			KRATOS_CHECK_EQUAL(triangle_splitter.mNegativeInterfacesParentIds[0], 0);

			// Check exterior faces
			KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 2);
			KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 3);

			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 0);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 0);

			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 1);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 1);

			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Y(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Y(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Y(), 0.0, tolerance);
		}

		KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3NoDivision, KratosCoreFastSuite)
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

			// Set the elemental distances vector
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.PointsNumber(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the triangle splitting utility
			DivideTriangle2D3 triangle_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			triangle_splitter.GenerateDivision();

			// Check general splitting values
			KRATOS_CHECK_IS_FALSE(triangle_splitter.mIsSplit);
			KRATOS_CHECK_EQUAL(triangle_splitter.mDivisionsNumber, 1);
			KRATOS_CHECK_EQUAL(triangle_splitter.mSplitEdgesNumber, 0);

		}
	}
}  // namespace Kratos.
