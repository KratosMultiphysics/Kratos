import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.unv_output_process as unv_output_process

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def ParseUnvBlocks(file_name):
    """Light UNV parser: returns a list of (dataset_id, [tokens]) blocks."""
    with open(file_name, "r") as f:
        tokens = f.read().split()

    blocks = []
    i = 0
    while i < len(tokens):
        if tokens[i] == "-1":
            i += 1
            if i >= len(tokens) or tokens[i] == "-1":
                continue
            dataset_id = int(tokens[i])
            i += 1
            content = []
            while i < len(tokens) and tokens[i] != "-1":
                content.append(tokens[i])
                i += 1
            i += 1  # skip closing "-1"
            blocks.append((dataset_id, content))
        else:
            i += 1
    return blocks

def SetUpModelPart(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
    model_part.SetBufferSize(1)

    # A small mesh with two linear triangles
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

    properties = model_part.CreateNewProperties(1)
    model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
    model_part.CreateNewElement("Element2D3N", 2, [1, 3, 4], properties)

    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, 25.0)

    sub_model_part = model_part.CreateSubModelPart("Skin")
    sub_model_part.AddNodes([1, 2, 3])
    sub_model_part.AddElements([1])

    model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.0
    model_part.ProcessInfo[KratosMultiphysics.STEP] = 1

class TestUnvOutputProcess(KratosUnittest.TestCase):
    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("Main.unv")

    def test_unv_output_process(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Main")
        SetUpModelPart(model_part)

        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"                    : "Main",
                "save_output_files_in_folder"        : false,
                "output_sub_model_parts"             : true,
                "nodal_solution_step_data_variables" : ["TEMPERATURE"]
            }
        }""")

        process = unv_output_process.Factory(settings, model)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.PrintOutput()

        self.assertTrue(os.path.isfile("Main.unv"))
        blocks = ParseUnvBlocks("Main.unv")
        dataset_ids = [block[0] for block in blocks]

        # Nodes (2411) and elements (2412) datasets must be present
        self.assertIn(2411, dataset_ids)
        self.assertIn(2412, dataset_ids)
        # Groups (2467) since output_sub_model_parts is true
        self.assertIn(2467, dataset_ids)
        # At least one nodal result dataset (2414) for TEMPERATURE
        self.assertIn(2414, dataset_ids)

        # 4 nodes -> 7 tokens each in the 2411 dataset
        nodes_block = next(content for did, content in blocks if did == 2411)
        self.assertEqual(len(nodes_block) // 7, 4)

if __name__ == "__main__":
    KratosUnittest.main()
