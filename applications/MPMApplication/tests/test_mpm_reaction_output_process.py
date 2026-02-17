import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.MPMApplication.mpm_grid_conforming_reaction_output_process as mpm_grid_conforming_reaction_output_process
import KratosMultiphysics.MPMApplication.mpm_non_conforming_reaction_output_process as mpm_non_conforming_reaction_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os
import pathlib

def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

class TestMPMReactionOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        super().setUp()

        # Create model
        self.model = KratosMultiphysics.Model()

        # Create model parts
        initial_mesh_model_part = self.model.CreateModelPart("InitialMesh")
        self.grid_model_part = self.model.CreateModelPart("Background_Grid")
        self.mpm_model_part = self.model.CreateModelPart("MPMModelPart")

        # Define Initial Mesh Model Part (used )
        initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        initial_mesh_sub_model_part = initial_mesh_model_part.CreateSubModelPart("sub_initial_mesh")

        # Define Background Grid model part
        self.grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        self.grid_model_part.CreateNewNode(0,  -0.5, -0.5, 0.0)
        self.grid_model_part.CreateNewNode(1,  -0.5, +0.5, 0.0)
        self.grid_model_part.CreateNewNode(2,  +0.5, -0.5, 0.0)
        self.grid_model_part.CreateNewNode(3,  +0.5, +0.5, 0.0)
        # Elements sub model part (needed for finding the background grid element containing a mp condition)
        grid_sub_model_part = self.grid_model_part.CreateSubModelPart("sub_background_grid")
        grid_sub_model_part.SetNodes(self.grid_model_part.GetNodes())
        grid_sub_model_part.CreateNewElement("Element2D4N", 1, [0, 1, 3, 2], grid_sub_model_part.GetProperties()[1])
        # Line conditions grid model part
        boundary_grid_model_part = self.grid_model_part.CreateSubModelPart("grid_conforming_reaction")
        boundary_grid_model_part.AddNodes([0, 1, 3])
        boundary_grid_model_part.CreateNewCondition("LineCondition2D2N", 1, [0, 1], self.grid_model_part.GetProperties()[1])
        boundary_grid_model_part.CreateNewCondition("LineCondition2D2N", 2, [1, 3], self.grid_model_part.GetProperties()[1])
        # Material point conditions mpm model part
        boundary_mpm_model_part = self.grid_model_part.CreateSubModelPart("non_conforming_reaction")
        boundary_mpm_model_part.AddNodes([0, 2, 3])
        boundary_mpm_model_part.CreateNewCondition("LineCondition2D2N", 3, [3, 2], self.grid_model_part.GetProperties()[1])
        boundary_mpm_model_part.CreateNewCondition("LineCondition2D2N", 4, [2, 0], self.grid_model_part.GetProperties()[1])
        for condition in boundary_mpm_model_part.Conditions:
            condition.Set(KratosMultiphysics.BOUNDARY, True)
            condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, 2)
            condition.SetValue(KratosMPM.MPC_BOUNDARY_CONDITION_TYPE, 1)

        # Define Material Point Model Part (this will contain mp elements and mp conditions)
        self.mpm_model_part.ProcessInfo = self.grid_model_part.ProcessInfo
        self.mpm_model_part.SetNodes(self.grid_model_part.GetNodes())
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)
        KratosMPM.GenerateMaterialPointElement(self.grid_model_part, initial_mesh_model_part, self.mpm_model_part, False)
        KratosMPM.GenerateMaterialPointCondition(self.grid_model_part, initial_mesh_model_part, self.mpm_model_part)

    def tearDown(self):
       kratos_utils.DeleteDirectoryIfExisting("test_mpm_reaction_output_process")

    def _set_solution(self):
        time = self.grid_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step = self.grid_model_part.ProcessInfo[KratosMultiphysics.STEP]

        for node in self.grid_model_part.Nodes:
            node_id = node.Id
            reaction_value = time + step + node_id
            reaction_value = [reaction_value*0.5, reaction_value*1.0, reaction_value*2.0]
            node.SetSolutionStepValue(KratosMultiphysics.REACTION, reaction_value)

        for i, cond in enumerate(self.mpm_model_part.GetSubModelPart("non_conforming_reaction").Conditions):
            reaction_value = time + step + i
            reaction_value = [reaction_value*0.25, reaction_value*0.5, reaction_value*1.0]
            cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_CONTACT_FORCE, [reaction_value], self.mpm_model_part.ProcessInfo)

    def test_mpm_grid_conforming_reaction_output_process(self):

        process_params = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"      : "Background_Grid.grid_conforming_reaction",
                "output_control_type"  : "step",
                "output_interval"      : 1,
                "print_format"         : ".8f",
                "output_file_settings" : {
                    "output_path" : "test_mpm_reaction_output_process",
                    "file_name" : "grid_conforming_reaction",
                    "file_extension" : "dat"
                }
            }
        }""")
        output_process = mpm_grid_conforming_reaction_output_process.Factory(process_params, self.model)

        time = 0.0
        dt = 0.2
        step = 0
        end_time = 0.6

        output_process.ExecuteInitialize()
        output_process.ExecuteBeforeSolutionLoop()
        while (time < end_time):
            time += dt
            self.mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            self.mpm_model_part.CloneTimeStep(time)
            self._set_solution()
            output_process.ExecuteInitializeSolutionStep()
            output_process.ExecuteFinalizeSolutionStep()
            if output_process.IsOutputStep():
                output_process.ExecuteBeforeOutputStep()
                output_process.PrintOutput()
                output_process.ExecuteAfterOutputStep()
        output_process.ExecuteFinalize()

        file_name  = process_params["Parameters"]["output_file_settings"]["file_name"].GetString()
        file_name += "." + process_params["Parameters"]["output_file_settings"]["file_extension"].GetString()
        output_path = pathlib.Path(process_params["Parameters"]["output_file_settings"]["output_path"].GetString())
        output_file = output_path/file_name
        reference_files_path = pathlib.Path("mpm_reaction_output_process_files")
        reference_file = reference_files_path/file_name
        params = KratosMultiphysics.Parameters("""{
           "reference_file_name" : "",
           "output_file_name"    : ""
        }""")
        params["reference_file_name"].SetString(str(GetFilePath(reference_file)))
        params["output_file_name"].SetString(str(output_file))
        CompareTwoFilesCheckProcess(params).Execute()

    def test_mpm_non_conforming_reaction_output_process(self):

        process_params = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"      : "MPMModelPart.non_conforming_reaction",
                "output_control_type"  : "step",
                "output_interval"      : 1,
                "print_format"         : ".8f",
                "output_file_settings" : {
                    "output_path" : "test_mpm_reaction_output_process",
                    "file_name" : "non_conforming_reaction",
                    "file_extension" : "dat"
                }
            }
        }""")
        output_process = mpm_non_conforming_reaction_output_process.Factory(process_params, self.model)

        time = 0.0
        dt = 0.2
        step = 0
        end_time = 0.6

        output_process.ExecuteInitialize()
        output_process.ExecuteBeforeSolutionLoop()
        while (time < end_time):
            time += dt
            self.mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            self.mpm_model_part.CloneTimeStep(time)
            self._set_solution()
            output_process.ExecuteInitializeSolutionStep()
            output_process.ExecuteFinalizeSolutionStep()
            if output_process.IsOutputStep():
                output_process.ExecuteBeforeOutputStep()
                output_process.PrintOutput()
                output_process.ExecuteAfterOutputStep()
        output_process.ExecuteFinalize()

        file_name  = process_params["Parameters"]["output_file_settings"]["file_name"].GetString()
        file_name += "." + process_params["Parameters"]["output_file_settings"]["file_extension"].GetString()
        output_path = pathlib.Path(process_params["Parameters"]["output_file_settings"]["output_path"].GetString())
        output_file = output_path/file_name
        reference_files_path = pathlib.Path("mpm_reaction_output_process_files")
        reference_file = reference_files_path/file_name
        params = KratosMultiphysics.Parameters("""{
           "reference_file_name" : "",
           "output_file_name"    : ""
        }""")
        params["reference_file_name"].SetString(str(GetFilePath(reference_file)))
        params["output_file_name"].SetString(str(output_file))
        CompareTwoFilesCheckProcess(params).Execute()

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.INFO)
    KratosUnittest.main()
