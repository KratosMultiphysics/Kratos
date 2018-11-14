//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  			     Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
// #include "includes/gid_io.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/embedded_skin_utility.h"

namespace Kratos
{
namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinUtility2D, KratosCoreFastSuite)
    {
        Model current_model;

        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
		Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
		Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.00, 10.00, 0.00);
		Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
		Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 10.00, 0.00, 0.00);

		Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   12,
			"element_name":     "Element2D3N"
		})");

		ModelPart &surface_part = current_model.CreateModelPart("Volume");
		surface_part.AddNodalSolutionStepVariable(DISTANCE);
		StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();

		// Generate the skin
		ModelPart &skin_part = current_model.CreateModelPart("Skin");
		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.4, 2.4, 0.0);
		skin_part.CreateNewNode(902, 7.6, 2.4, 0.0);
		skin_part.CreateNewNode(903, 7.6, 7.6, 0.0);
		skin_part.CreateNewNode(904, 2.4, 7.6, 0.0);
		skin_part.CreateNewNode(905, 3.9, 3.9, 0.0);
		skin_part.CreateNewNode(906, 6.1, 3.9, 0.0);
		skin_part.CreateNewNode(907, 6.1, 6.1, 0.0);
		skin_part.CreateNewNode(908, 3.9, 6.1, 0.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element2D2N", 901, {{901,902}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 902, {{902,903}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 903, {{903,904}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 904, {{904,901}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 905, {{905,906}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 906, {{906,907}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 907, {{907,908}}, p_properties);
		skin_part.CreateNewElement("Element2D2N", 908, {{908,905}}, p_properties);

		// Compute distance
		CalculateDistanceToSkinProcess<2>(surface_part, skin_part).Execute();

        // Generate the skin
        ModelPart &generated_skin_part = current_model.CreateModelPart("GeneratedSkinPart");
        EmbeddedSkinUtility<2> embedded_skin_utility(surface_part, generated_skin_part, "continuous");
		embedded_skin_utility.GenerateSkin();

        KRATOS_CHECK_EQUAL(generated_skin_part.NumberOfNodes(), 152);
        KRATOS_CHECK_EQUAL(generated_skin_part.NumberOfConditions(), 76);

        // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/surface_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		// gid_io_fluid.InitializeMesh(0.00);
		// gid_io_fluid.WriteMesh(surface_part.GetMesh());
		// gid_io_fluid.FinalizeMesh();
		// gid_io_fluid.InitializeResults(0, surface_part.GetMesh());
		// gid_io_fluid.WriteNodalResults(DISTANCE, surface_part.Nodes(), 0, 0);
		// gid_io_fluid.FinalizeResults();

		// GidIO<> gid_io_skin("/home/rzorrilla/Desktop/skin_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		// gid_io_skin.InitializeMesh(0.00);
		// gid_io_skin.WriteMesh(skin_part.GetMesh());
		// gid_io_skin.FinalizeMesh();
		// gid_io_skin.InitializeResults(0, skin_part.GetMesh());
		// gid_io_skin.FinalizeResults();

		// GidIO<> gid_io_generated_skin("/home/rzorrilla/Desktop/generated_skin_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
		// gid_io_generated_skin.InitializeMesh(0.00);
		// gid_io_generated_skin.WriteMesh(generated_skin_part.GetMesh());
		// gid_io_generated_skin.FinalizeMesh();
		// gid_io_generated_skin.InitializeResults(0, generated_skin_part.GetMesh());
		// gid_io_generated_skin.FinalizeResults(); 
    }

}
}  // namespace Kratos.
