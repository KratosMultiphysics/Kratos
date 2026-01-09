import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.MPMApplication.mpm_write_energy_output_process as mpm_write_energy_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os
import pathlib

def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

class TestMPMWriteEnergyOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        super().setUp()

        # Create model
        self.model = KratosMultiphysics.Model()

        # Create model parts
        initial_mesh_model_part = self.model.CreateModelPart("InitialMesh")
        grid_model_part = self.model.CreateModelPart("Background_Grid")
        mpm_model_part = self.model.CreateModelPart("MPMModelPart")

        # Define Initial Mesh Model Part
        initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        initial_mesh_sub_model_part = initial_mesh_model_part.CreateSubModelPart("SubInitialMesh")
        initial_mesh_sub_model_part.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
        initial_mesh_sub_model_part.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 2)
        initial_mesh_sub_model_part.CreateNewNode(1, -0.5, -0.5, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(2, -0.5,  0.0, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(3, -0.5, +0.5, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(4,  0.0, -0.5, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(5,  0.0,  0.0, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(6,  0.0, +0.5, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(7, +0.5, -0.5, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(8, +0.5,  0.0, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(9, +0.5, +0.5, 0.0)
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 1, [1, 4, 5, 2], initial_mesh_sub_model_part.GetProperties()[1])
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 2, [2, 5, 6, 3], initial_mesh_sub_model_part.GetProperties()[1])
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 3, [4, 7, 8, 5], initial_mesh_sub_model_part.GetProperties()[1])
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 4, [5, 8, 9, 6], initial_mesh_sub_model_part.GetProperties()[1])

        # Define Background Grid model part
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        grid_sub_model_part = grid_model_part.CreateSubModelPart("SubBackgroundGrid")
        grid_sub_model_part.CreateNewNode(1,  -0.5, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(2,  -0.5, +0.5, 0.0)
        grid_sub_model_part.CreateNewNode(4,  +0.5, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(5,  +0.5, +0.5, 0.0)
        grid_sub_model_part.CreateNewElement("Element2D4N", 1, [ 1,  2,  5,  4], grid_sub_model_part.GetProperties()[1])

        # Define Material Point Model Part (this will contain mp elements)
        mpm_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        mpm_model_part.SetNodes(grid_model_part.GetNodes())
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)
        KratosMPM.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, mpm_model_part, False)
        KratosMPM.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, mpm_model_part)

        for elem in mpm_model_part.Elements:
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY, [[0.1, 0.2, 0.3]], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION, [[0.5, 0.5, 0.5]], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_CAUCHY_STRESS_VECTOR, [KratosMultiphysics.Vector([1, 2, 3, 4, 5, 6])], 1, mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_ALMANSI_STRAIN_VECTOR, [KratosMultiphysics.Vector([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])], 1, mpm_model_part.ProcessInfo)

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_mpm_write_energy")

    def set_solution(self):
        mpm_model_part = self.model.GetModelPart("MPMModelPart")
        time = mpm_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step = mpm_model_part.ProcessInfo[KratosMultiphysics.STEP]

        for i, elem in enumerate(mpm_model_part.Elements):
            coord = elem.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, mpm_model_part.ProcessInfo)[0]
            velocity = elem.CalculateOnIntegrationPoints(KratosMPM.MP_VELOCITY, mpm_model_part.ProcessInfo)[0]
            volume_acceleration = elem.CalculateOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION, mpm_model_part.ProcessInfo)[0]
            mass = elem.CalculateOnIntegrationPoints(KratosMPM.MP_MASS, mpm_model_part.ProcessInfo)[0]
            volume = elem.CalculateOnIntegrationPoints(KratosMPM.MP_VOLUME, mpm_model_part.ProcessInfo)[0]
            cauchy_stress = elem.CalculateOnIntegrationPoints(KratosMPM.MP_CAUCHY_STRESS_VECTOR, mpm_model_part.ProcessInfo)[0]
            almansi_strain = elem.CalculateOnIntegrationPoints(KratosMPM.MP_ALMANSI_STRAIN_VECTOR, mpm_model_part.ProcessInfo)[0]
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_COORD, [[coord[0]*time, coord[1]+i*step*0.1, 0.1*i+coord[2]]], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY, [[velocity[0]+step, velocity[1]+i*0.2, velocity[2]]], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION, [[volume_acceleration[0]-0.1, volume_acceleration[1]+i*0.05, 1.1*volume_acceleration[2]]], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_MASS, [mass*0.3*(i+1)], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME, [volume+0.2*i], mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_CAUCHY_STRESS_VECTOR, [(time+1)*KratosMultiphysics.Vector([1, 2, 3, 4, 5, 6])], 6, mpm_model_part.ProcessInfo)
            elem.SetValuesOnIntegrationPoints(KratosMPM.MP_ALMANSI_STRAIN_VECTOR, [(1+0.1*step)*KratosMultiphysics.Vector([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])], 6, mpm_model_part.ProcessInfo)

    def test_mpm_write_energy_output_process(self):

        mpm_write_energy_output_parameters = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"      : "MPMModelPart",
                "interval"             : [0.0, 1e30],
                "print_format"         : ".8f",
                "output_file_settings" : {
                    "output_path" : "test_mpm_write_energy",
                    "file_extension" : "dat"
                }
            }
        }""")
        output_process = mpm_write_energy_output_process.Factory(mpm_write_energy_output_parameters, self.model)

        mpm_model_part = self.model.GetModelPart("MPMModelPart")
        time = 0.0
        dt = 0.2
        step = 0
        end_time = 0.6

        output_process.ExecuteInitialize()
        output_process.ExecuteBeforeSolutionLoop()
        while (time < end_time):
            time += dt
            mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            mpm_model_part.CloneTimeStep(time)
            self.set_solution()
            output_process.ExecuteInitializeSolutionStep()
            output_process.ExecuteFinalizeSolutionStep()
            if output_process.IsOutputStep():
                output_process.ExecuteBeforeOutputStep()
                output_process.PrintOutput()
                output_process.ExecuteAfterOutputStep()
        output_process.ExecuteFinalize()

        file_name = "MPMModelPart_energy.dat"
        output_path = pathlib.Path(mpm_write_energy_output_parameters["Parameters"]["output_file_settings"]["output_path"].GetString())
        output_file = output_path/file_name
        reference_files_path = pathlib.Path("mpm_write_energy_output_process_file")
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
