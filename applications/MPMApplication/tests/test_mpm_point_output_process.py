import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.MPMApplication.mpm_multiple_points_output_process as mpm_multiple_points_output_process
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import pathlib

def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

class TestMPMPointOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        super().setUp()

        self.model = KratosMultiphysics.Model()
        initial_mesh_model_part = self.model.CreateModelPart("InitialMesh")
        self.grid_model_part = self.model.CreateModelPart("Background_Grid")
        self.mpm_model_part = self.model.CreateModelPart("MPMModelPart")

        # Define Initial Mesh Model Part (used for defining material point elements)
        initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        initial_mesh_sub_model_part = initial_mesh_model_part.CreateSubModelPart("SubInitialMesh")
        # Number of material point elements for each element of the initial mesh
        initial_mesh_sub_model_part.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
        # Nodes for the initial mesh
        initial_mesh_sub_model_part.CreateNewNode(1, -0.25, -0.25, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(2, -0.25,  0.00, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(3, -0.25, +0.25, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(4,  0.00, -0.25, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(5,  0.00,  0.00, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(6,  0.00, +0.25, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(7, +0.25, -0.25, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(8, +0.25,  0.00, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(9, +0.25, +0.25, 0.0)
        # Elements for initial mesh
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 5, 4], initial_mesh_sub_model_part.GetProperties()[1])
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 6, 5], initial_mesh_sub_model_part.GetProperties()[1])
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 3, [4, 5, 8, 7], initial_mesh_sub_model_part.GetProperties()[1])
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 9, 8], initial_mesh_sub_model_part.GetProperties()[1])

        # Define Background Grid model part
        self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        grid_sub_model_part = self.grid_model_part.CreateSubModelPart("SubBackgroundGrid")
        # Nodes for background grid
        grid_sub_model_part.CreateNewNode(1,  -1.0, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(2,  -1.0,  0.0, 0.0)
        grid_sub_model_part.CreateNewNode(3,  -1.0, +0.5, 0.0)
        grid_sub_model_part.CreateNewNode(4,  -0.5, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(5,  -0.5,  0.0, 0.0)
        grid_sub_model_part.CreateNewNode(6,  -0.5, +0.5, 0.0)
        grid_sub_model_part.CreateNewNode(7,   0.0, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(8,   0.0,  0.0, 0.0)
        grid_sub_model_part.CreateNewNode(9,   0.0, +0.5, 0.0)
        grid_sub_model_part.CreateNewNode(10, +0.5, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(11, +0.5,  0.0, 0.0)
        grid_sub_model_part.CreateNewNode(12, +0.5, +0.5, 0.0)
        grid_sub_model_part.CreateNewNode(13, +1.0, -0.5, 0.0)
        grid_sub_model_part.CreateNewNode(14, +1.0,  0.0, 0.0)
        grid_sub_model_part.CreateNewNode(15, +1.0, +0.5, 0.0)
        # Elements for background grid
        grid_sub_model_part.CreateNewElement("Element2D4N", 1, [ 1,  2,  5,  4], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 2, [ 2,  3,  6,  5], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 3, [ 4,  5,  8,  7], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 4, [ 5,  6,  9,  8], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 5, [ 7,  8, 11, 10], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 6, [ 8,  9, 12, 11], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 7, [10, 11, 14, 13], grid_sub_model_part.GetProperties()[1])
        grid_sub_model_part.CreateNewElement("Element2D4N", 8, [11, 12, 15, 14], grid_sub_model_part.GetProperties()[1])

        # Interface
        grid_interface = self.grid_model_part.CreateSubModelPart("InterfaceConditions")
        # Nodes for interface condition
        grid_interface.CreateNewNode(16, -0.2, 0.0, 0.0)
        grid_interface.CreateNewNode(17, -0.1, 0.0, 0.0)
        grid_interface.CreateNewNode(18,  0.0, 0.0, 0.0)
        grid_interface.CreateNewNode(19, +0.1, 0.0, 0.0)
        grid_interface.CreateNewNode(20, +0.2, 0.0, 0.0)
        # Conditions
        grid_interface.CreateNewCondition("LineCondition2D2N", 1, [16, 17], grid_interface.GetProperties()[1])
        grid_interface.CreateNewCondition("LineCondition2D2N", 2, [17, 18], grid_interface.GetProperties()[1])
        grid_interface.CreateNewCondition("LineCondition2D2N", 3, [18, 19], grid_interface.GetProperties()[1])
        grid_interface.CreateNewCondition("LineCondition2D2N", 4, [19, 20], grid_interface.GetProperties()[1])
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, grid_interface.Conditions)
        for condition in grid_interface.Conditions:
            condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, 1)
            condition.SetValue(KratosMPM.MPC_BOUNDARY_CONDITION_TYPE, 1)

        # Define Material Point Model Part (this will contain mp elements)
        self.mpm_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Activate elements of initial mesh model part
        self.mpm_model_part.SetNodes(self.grid_model_part.GetNodes())
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)
        # Generate Material Point Elements and Conditions
        KratosMPM.GenerateMaterialPointElement(self.grid_model_part, initial_mesh_model_part, self.mpm_model_part, False)
        KratosMPM.GenerateMaterialPointCondition(self.grid_model_part, initial_mesh_model_part, self.mpm_model_part)

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("test_material_point_output")

    def _set_solution(self):
        t = self.mpm_model_part.ProcessInfo[KratosMultiphysics.TIME] + 0.150
        s = self.mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] + 0.2

        for element in self.mpm_model_part.Elements:
            coord = element.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, self.mpm_model_part.ProcessInfo)[0]
            updated_coord = [coord[0]*1.25, coord[1]+1.75, coord[2]]
            element.SetValuesOnIntegrationPoints(KratosMPM.MP_COORD, [updated_coord], self.mpm_model_part.ProcessInfo)

            density = element.CalculateOnIntegrationPoints(KratosMPM.MP_DENSITY, self.mpm_model_part.ProcessInfo)[0]
            updated_density = density + 0.2
            element.SetValuesOnIntegrationPoints(KratosMPM.MP_DENSITY, [updated_density], self.mpm_model_part.ProcessInfo)

            displacement = element.CalculateOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, self.mpm_model_part.ProcessInfo)[0]
            updated_displacement = [displacement[0]+1.5*s, displacement[1]+s/2, displacement[2]]
            element.SetValuesOnIntegrationPoints(KratosMPM.MP_DISPLACEMENT, [updated_displacement], self.mpm_model_part.ProcessInfo)

        for condition in self.mpm_model_part.Conditions:
            coord = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, self.mpm_model_part.ProcessInfo)[0]
            updated_coord = [coord[0]+0.35, coord[1]+1+t, coord[2]]
            condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_COORD, [updated_coord], self.mpm_model_part.ProcessInfo)

            displacement = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, self.mpm_model_part.ProcessInfo)[0]
            updated_displacement = [displacement[0]+t, 1.1*displacement[1]-s, displacement[2]]
            condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, [updated_displacement], self.mpm_model_part.ProcessInfo)

            area = condition.CalculateOnIntegrationPoints(KratosMPM.MPC_AREA, self.mpm_model_part.ProcessInfo)[0]
            updated_area = (area+t) * 1.75
            condition.SetValuesOnIntegrationPoints(KratosMPM.MPC_AREA, [updated_area], self.mpm_model_part.ProcessInfo)

    def _get_process_parameters(self, entity_type):
        parameters = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"      : "MPMModelPart",
                "entity_type"          : "",
                "output_control_type"  : "step",
                "output_interval"      : 1.0,
                "positions"            : [],
                "output_variables"     : [],
                "search_tolerance"     : 1e-3,
                "print_format"         : ".5e",
                "output_file_settings" : {
                    "file_name"         : "",
                    "output_path"       : "test_material_point_output",
                    "file_extension"    : "dat"
                }
            }
        }""")
        parameters["Parameters"]["output_file_settings"]["file_name"].SetString(f"mp_process_output_{entity_type}")
        parameters["Parameters"]["entity_type"].SetString(entity_type)
        return parameters

    def _run_and_check_output_process(self, process_params):
        """Run the output process through a time loop and compare results to reference"""
        test_process = mpm_multiple_points_output_process.Factory(process_params, self.model)

        time = 0.0
        dt = 0.2
        end_time = 1.0

        test_process.ExecuteInitialize()
        test_process.ExecuteBeforeSolutionLoop()
        while (time < end_time):
            time += dt
            self.mpm_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            self.mpm_model_part.CloneTimeStep(time)
            self._set_solution()
            test_process.ExecuteInitializeSolutionStep()
            test_process.ExecuteFinalizeSolutionStep()
            if test_process.IsOutputStep():
                test_process.ExecuteBeforeOutputStep()
                test_process.PrintOutput()
                test_process.ExecuteAfterOutputStep()
        test_process.ExecuteFinalize()

        # Compare output file with reference file
        reference_files_path = pathlib.Path("material_point_output_process_files")
        file_name = process_params["Parameters"]["output_file_settings"]["file_name"].GetString()
        output_path = process_params["Parameters"]["output_file_settings"]["output_path"].GetString()
        file_extension = process_params["Parameters"]["output_file_settings"]["file_extension"].GetString()
        self._check(f"{file_name}_1.{file_extension}", output_path, reference_files_path)
        self._check(f"{file_name}_2.{file_extension}", output_path, reference_files_path)

    def _check(self, file_name, output_path, reference_files_path):
        output_file = pathlib.Path(output_path)/pathlib.Path(file_name)
        reference_file = GetFilePath(reference_files_path / file_name)
        params = KratosMultiphysics.Parameters()
        params.AddString("reference_file_name", str(reference_file))
        params.AddString("output_file_name", str(output_file))
        params.AddString("comparison_type", "deterministic")
        CompareTwoFilesCheckProcess(params).Execute()

    def test_mpm_output_process_element_2D(self):
        parameters = self._get_process_parameters("element")
        positions = KratosMultiphysics.Matrix([[-0.125,-0.125,0.0],[-0.12499,0.12501,0.0]])
        parameters["Parameters"]["positions"].SetMatrix(positions)
        parameters["Parameters"]["output_variables"].SetStringArray(["MP_COORD","MP_DISPLACEMENT","MP_DENSITY"])
        self._run_and_check_output_process(parameters)

    def test_mpm_output_process_condition_2D(self):
        parameters = self._get_process_parameters("condition")
        positions = KratosMultiphysics.Matrix([[-0.0501, 0.0, 0.0], [0.15, -0.00001, 0.0]])
        parameters["Parameters"]["positions"].SetMatrix(positions)
        parameters["Parameters"]["output_variables"].SetStringArray(["MPC_COORD","MPC_DISPLACEMENT","MPC_AREA"])
        self._run_and_check_output_process(parameters)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.DETAIL)
    KratosUnittest.main()
