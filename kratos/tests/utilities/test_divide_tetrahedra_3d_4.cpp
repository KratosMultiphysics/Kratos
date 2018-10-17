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

namespace Kratos
{
	namespace Testing
	{

		KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4Horizontal, KratosCoreFastSuite)
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
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 4> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the tetrahedra splitting utility
			DivideTetrahedra3D4 tetrahedra_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			tetrahedra_splitter.GenerateDivision();

			// Call the intersection generation method
			tetrahedra_splitter.GenerateIntersectionsSkin();

			// Call the positive exterior faces generation method
			std::vector < unsigned int > pos_ext_faces_parent_ids;
			std::vector < DivideTetrahedra3D4::IndexedPointGeometryPointerType > pos_ext_faces;
			tetrahedra_splitter.GenerateExteriorFaces(
				pos_ext_faces,
				pos_ext_faces_parent_ids,
				tetrahedra_splitter.mPositiveSubdivisions);

			// Call the negative exterior faces generation method
			std::vector < unsigned int > neg_ext_faces_parent_ids;
			std::vector < DivideTetrahedra3D4::IndexedPointGeometryPointerType > neg_ext_faces;
			tetrahedra_splitter.GenerateExteriorFaces(
				neg_ext_faces,
				neg_ext_faces_parent_ids,
				tetrahedra_splitter.mNegativeSubdivisions);

			const double tolerance = 1e-10;

			// Check general splitting values
			KRATOS_CHECK(tetrahedra_splitter.mIsSplit);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 4);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 3);

			// Check split edges
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[0],  0);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[1],  1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[2],  2);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[3],  3);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[4], -1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[5], -1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[6],  6);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[7], -1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[8],  8);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[9],  9);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[10],-1);

			// Check subdivisions
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[2].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[3].Z(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[2].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[3].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[0].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[2].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[3].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[3].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[2].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[3].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[3].Z(), 0.5, tolerance);

			// Check interfaces
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mPositiveInterfacesParentIds[0], 0);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mNegativeInterfacesParentIds[0], 2);

			// Check exterior faces
			KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 3);
			KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 7);

			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 0);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 0);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[2], 0);

			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 1);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 2);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[3], 2);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[4], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[5], 1);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[6], 0);

			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Z(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Z(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].Z(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[0].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[2].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[6])[2].Z(), 0.0, tolerance);
		}
		
		KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4Oblique, KratosCoreFastSuite)
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
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 4> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);
			
			// Build the tetrahedra splitting utility
			DivideTetrahedra3D4 tetrahedra_splitter(r_geometry, r_elemental_distances);
			
			// Call the divide geometry method
			tetrahedra_splitter.GenerateDivision();

			// Call the intersection generation method
			tetrahedra_splitter.GenerateIntersectionsSkin();

			// Call the positive exterior faces generation method
			std::vector < unsigned int > pos_ext_faces_parent_ids;
			std::vector < DivideTetrahedra3D4::IndexedPointGeometryPointerType > pos_ext_faces;
			tetrahedra_splitter.GenerateExteriorFaces(
				pos_ext_faces,
				pos_ext_faces_parent_ids,
				tetrahedra_splitter.mPositiveSubdivisions);

			// Call the negative exterior faces generation method
			std::vector < unsigned int > neg_ext_faces_parent_ids;
			std::vector < DivideTetrahedra3D4::IndexedPointGeometryPointerType > neg_ext_faces;
			tetrahedra_splitter.GenerateExteriorFaces(
				neg_ext_faces,
				neg_ext_faces_parent_ids,
				tetrahedra_splitter.mNegativeSubdivisions);

			const double tolerance = 1e-10;

			// Check general splitting values
			KRATOS_CHECK(tetrahedra_splitter.mIsSplit);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 6);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 4);
			
			// Check split edges
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[0],  0);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[1],  1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[2],  2);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[3],  3);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[4],  4);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[5], -1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[6],  6);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[7],  7);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[8], -1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[9],  9);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdges[10],-1);
			
			// Check subdivisions
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[2].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[0])[3].Z(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[2].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[1])[3].Z(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[2].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveSubdivisions[2])[3].Z(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[2].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[3].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[0])[3].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[2].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[3].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[1])[3].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[0].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[2].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[3].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[3].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeSubdivisions[2])[3].Z(), 0.5, tolerance);

			// Check interfaces
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[0])[2].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mPositiveInterfaces[1])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[0])[2].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*tetrahedra_splitter.mNegativeInterfaces[1])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mPositiveInterfacesParentIds[0], 0);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mPositiveInterfacesParentIds[1], 1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mNegativeInterfacesParentIds[0], 1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mNegativeInterfacesParentIds[1], 2);

			// Check exterior faces
			KRATOS_CHECK_EQUAL(pos_ext_faces.size(), 6);
			KRATOS_CHECK_EQUAL(neg_ext_faces.size(), 6);

			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[0], 1);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[1], 2);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[2], 0);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[3], 0);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[4], 2);
			KRATOS_CHECK_EQUAL(pos_ext_faces_parent_ids[5], 2);

			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[0], 2);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[1], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[2], 1);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[3], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[4], 0);
			KRATOS_CHECK_EQUAL(neg_ext_faces_parent_ids[5], 2);

			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[1].Z(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[0])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[1].Z(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[1])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[1].Z(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[2])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[1].Z(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[3])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[1].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[4])[2].Z(), 1.0, tolerance);

			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[0].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[2].X(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*pos_ext_faces[5])[2].Z(), 0.0, tolerance);
			
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[0])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[1])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[0].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[1].Z(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[2])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[1].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[2].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[3])[2].Z(), 0.5, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[0].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[1].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[1].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[4])[2].Z(), 0.0, tolerance);

			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[0].X(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[0].Y(), 1.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[0].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[1].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[1].Y(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[1].Z(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[2].X(), 0.5, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[2].Y(), 0.0, tolerance);
			KRATOS_CHECK_NEAR((*neg_ext_faces[5])[2].Z(), 0.0, tolerance);
		}
		
		KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTetrahedra3D4NoDivision, KratosCoreFastSuite)
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
			base_model_part.Nodes()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
			base_model_part.Nodes()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
			base_model_part.Nodes()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;
			base_model_part.Nodes()[4].FastGetSolutionStepValue(DISTANCE) = 1.0;
			
			// Set the elemental distances vector
			Geometry < Node < 3 > >& r_geometry = base_model_part.Elements()[1].GetGeometry();

			array_1d<double, 4> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}

			base_model_part.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			Vector& r_elemental_distances = base_model_part.Elements()[1].GetValue(ELEMENTAL_DISTANCES);

			// Build the tetrahedra splitting utility
			DivideTetrahedra3D4 tetrahedra_splitter(r_geometry, r_elemental_distances);

			// Call the divide geometry method
			tetrahedra_splitter.GenerateDivision();

			// Check general splitting values
			KRATOS_CHECK_IS_FALSE(tetrahedra_splitter.mIsSplit);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mDivisionsNumber, 1);
			KRATOS_CHECK_EQUAL(tetrahedra_splitter.mSplitEdgesNumber, 0);

		}
	}
}  // namespace Kratos.
