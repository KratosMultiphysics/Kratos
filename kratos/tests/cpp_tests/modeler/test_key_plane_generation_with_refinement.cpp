//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_application.h"
#include "includes/kernel.h"

#include "modeler/voxel_mesh_generator_modeler.h"
#include "modeler/key_plane_generation/key_plane_generation_factory.h"

namespace Kratos::Testing {

namespace {
void WriteCubeSkinMeshMdpaFile()
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

        9    -0.03    -0.02    0.05
        10   -0.01    -0.02    0.07
        11   -0.01    -0.01    0.05
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
        13    1    9    10   11
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
Begin SubModelPart refinement1
Begin SubModelPartNodes
9
10
11
End SubModelPartNodes
Begin SubModelPartElements
13
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
    ModelPartIO model_part_io_write(output_file_name, IO::WRITE);
    model_part_io_write.WriteModelPart(model_part_0);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorModelerWithRefinement, KratosCoreFastSuite)
{
    using namespace Kratos;

	WriteCubeSkinMeshMdpaFile();

    Parameters mesher_parameters(R"(
    {
        "output_model_part_name" : "output_model_part",
        "input_model_part_name" : "skin_model_part",
        "mdpa_file_name" : "cube_skin_mesh",
        "key_plane_generator": {
            "type" : "outer_shell_with_refinement",
            "Parameters" : {
                "voxel_sizes" : [0.005, 0.005, 0.005],
                "refinement_zones":[
                {
                    "refined_modelpart" : "skin_model_part.refinement1",
                    "refined_zone" : {},
                    "voxel_sizes_ratio": [0.5, 0.5, 0.5]
                }]

            }
        }
    })");

    Model current_model;
    current_model.CreateModelPart("main_model_part");

    // Generate the skin
    current_model.CreateModelPart("skin_model_part");

    // Generating the mesh
    auto voxel_mesher = VoxelMeshGeneratorModeler(current_model, mesher_parameters);
    voxel_mesher.SetupGeometryModel();
    voxel_mesher.ReadModelParts();

    KeyPlaneGenerationFactory factory;
    auto p_key_plane_generator = factory.Create(voxel_mesher, mesher_parameters["key_plane_generator"]);
    p_key_plane_generator->ValidateParameters();
    p_key_plane_generator->Generate();

    std::vector<double> x_key_planes {-0.0325,-0.03,-0.0275,-0.025,-0.0225,-0.02,-0.0175,-0.015,-0.0125,-0.01,-0.0075,-0.005,-0.0025,0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325};
    std::vector<double> y_key_planes {-0.035,-0.03,-0.025,-0.02,-0.0175,-0.015,-0.0125,-0.01,-0.0075,-0.005,-0.0025,0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325};
    std::vector<double> z_key_planes {0.025,0.03,0.034,0.038,0.042,0.046,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.0725,0.075,0.0775,0.08,0.0825,0.085,0.0875,0.09,0.0925};

    // delete work file
    std::remove("cube_skin_mesh.mdpa");

    KRATOS_EXPECT_VECTOR_NEAR(x_key_planes, voxel_mesher.GetKeyPlanes(0), 1e-6);
    KRATOS_EXPECT_VECTOR_NEAR(y_key_planes, voxel_mesher.GetKeyPlanes(1), 1e-6);
    KRATOS_EXPECT_VECTOR_NEAR(z_key_planes, voxel_mesher.GetKeyPlanes(2), 1e-6);
}

} // namespace Kratos::Testing
