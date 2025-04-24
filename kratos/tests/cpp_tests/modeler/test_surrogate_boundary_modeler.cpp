//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes Danes
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_application.h"
#include "includes/kernel.h"

#include "modeler/surrogate_boundary_modeler.h"

namespace Kratos::Testing {

namespace {
void WriteCubeSkinMeshMdpaFileForSBModelerTest()
{
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
	Begin ModelPartData
End ModelPartData
Begin Properties 0
End Properties
Begin Properties 1
End Properties
Begin Nodes
		1     0.03     0.03    0.09
		2     0.03    -0.03    0.09
		3     0.03    -0.03    0.03
		4     0.03     0.03    0.03
		5    -0.03     0.03    0.09
		6    -0.03     0.03    0.03
		7    -0.03    -0.03    0.09
		8    -0.03    -0.03    0.03
End Nodes
Begin Elements Element3D3N
		1     1    1    2    3
		2     1    3    4    1
		3     1    5    1    4
		4     1    4    6    5
		5     1    7    5    6
		6     1    6    8    7
		7     1    7    8    3
		8     1    3    2    7
		9     1    7    2    1
		10    1    1    5    7
		11    1    8    6    4
		12    1    4    3    8
End Elements
Begin SubModelPart workpiece
Begin SubModelPartNodes
		1
		2
		3
		4
		5
		6
		7
		8
End SubModelPartNodes
Begin SubModelPartElements
		1
		2
		3
		4
		5
		6
		7
		8
		9
		10
		11
		12
End SubModelPartElements
End SubModelPart
	)input"));

    Model current_model;

    ModelPartIO model_part_io(p_input);

    ModelPart& model_part_0 =
        current_model.CreateModelPart("tmp_main_model_part");
    model_part_io.ReadModelPart(model_part_0);

    // Create the output .mdpa file
    std::string output_file_name = "cube_skin_mesh";
    std::fstream output_file;
    output_file.open(output_file_name + ".mdpa", std::fstream::out);
    output_file.close();

    // Fill the output .mdpa file
    ModelPartIO* model_part_io_write = new ModelPartIO(output_file_name, IO::WRITE);
    model_part_io_write->WriteModelPart(model_part_0);
}
} // namespace

	KRATOS_TEST_CASE_IN_SUITE(VoxelMesherWithVectorDistances, KratosCoreFastSuite)
	{
		using namespace Kratos;

		WriteCubeSkinMeshMdpaFileForSBModelerTest();
		std::cout << "Cube writen" << std::endl;

		Parameters mesher_parameters(R"(
			{
				"output_model_part_name" : "main_model_part",
				"input_model_part_name" : "skin_model_part",
				"mdpa_file_name" : "cube_skin_mesh",
				"key_plane_generator": {
					"Parameters" : {
						"voxel_sizes" : [0.025, 0.025, 0.025],
						"min_point" : [-0.05, -0.05, 0],
						"max_point" : [0.05, 0.05, 0.1]
					}
				},
				"coloring_settings_list": [
				{
					"type" : "cells_in_touch",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				},
				{
					"type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				}
				],
				"entities_generator_list": [
				{
					"type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
					"properties_id": 1
				} 
				]
			})");
			
		Model current_model;
		current_model.CreateModelPart("main_model_part");

		// Generate the skin
		current_model.CreateModelPart("skin_model_part");

		std::cout << "Modelpart created" << std::endl;

		// Generating the mesh
		auto voxel_mesher = SurrogateBoundaryModeler(current_model, mesher_parameters);
		voxel_mesher.SetupGeometryModel();
		voxel_mesher.SetupModelPart();

		std::cout << "SurrogateBoundaryModeler created" << std::endl;

		voxel_mesher.ComputeSurrogateBoundary();

		std::cout << "Distances and colors computed" << std::endl;
		
		auto& nodalSBdata = voxel_mesher.GetSurrogateBoundaryNodes();

		for (SurrogateBoundaryModeler::SurrogateBoundaryNode& node : nodalSBdata) 
		{
			// We could also use GetSurrogateBoundaryNode(i,j,k)
			if (node.IsActive()) 
			{
				std::cout << *node.GetNodePtr() 
						  << " \n  Vector distance to skin: " << node.GetVectorDistance()
						  << " \n  Signed distance to skin: " << node.GetSignedDistance() 
						  << " \n  Is inside: " << node.IsInside() << std::endl;
			
                if (std::abs(node.GetNodePtr()->X()) > 0.03 ||
                    std::abs(node.GetNodePtr()->Y()) > 0.03 ||
                    node.GetNodePtr()->Z() < 0.03 || node.GetNodePtr()->Z() > 0.09)
                {
                    KRATOS_CHECK_IS_FALSE(node.IsInside());
                }
            }
		}
	}

} // namespace Kratos::Testing
