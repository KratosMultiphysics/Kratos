import KratosMultiphysics

import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing analysis
from KratosMultiphysics.FreeSurfaceApplication.free_surface_analysis import FreeSurfaceAnalysis

# Importing post-process
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

class KratosFreeSurfaceGeneralTests(KratosUnittest.TestCase):
    def setUp(self):
        self.print_output = False
        self.print_results = False

    def test_reservoir_2d(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "reservoir_2d/reservoir_2d_results.json"
            parameters_filename = "reservoir_2d/reservoir_2d_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            model = KratosMultiphysics.Model()
            simulation = FreeSurfaceAnalysis(model, parameters)
            simulation.Run()

            # self._check_results(model_part, A, b)
            if self.print_results:
                self.__print_results(model, results_filename, parameters)
            if self.print_output:
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()))
            self.__check_results(model, results_filename, parameters)

    def test_reservoir_3d(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "reservoir_3d/reservoir_3d_results.json"
            parameters_filename = "reservoir_3d/reservoir_3d_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            model = KratosMultiphysics.Model()
            simulation = FreeSurfaceAnalysis(model, parameters)
            simulation.Run()

            # self._check_results(model_part, A, b)
            if self.print_results:
                self.__print_results(model, results_filename, parameters)
            if self.print_output:
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()))
            self.__check_results(model, results_filename, parameters)

    def __print_results(self, model, results_filename, parameters):
        json_output_settings = KratosMultiphysics.Parameters(r"""
        {
            "output_variables": ["DISTANCE"],
            "output_file_name": "",
            "model_part_name": "",
            "time_frequency": 0.0
        }""")
        json_output_settings["model_part_name"].SetString(parameters["solver_settings"]["model_part_name"].GetString() + '.Results')
        json_output_settings["output_file_name"].SetString(results_filename)
        self.json_output = JsonOutputProcess(model, json_output_settings)
        self.json_output.ExecuteInitialize()
        self.json_output.ExecuteBeforeSolutionLoop()
        self.json_output.ExecuteFinalizeSolutionStep()

    def __check_results(self, model, results_filename, parameters):
        from_json_check_result_settings = KratosMultiphysics.Parameters(r"""
        {
            "check_variables": ["DISTANCE"],
            "input_file_name": "",
            "model_part_name": ""
        }""")
        from_json_check_result_settings["model_part_name"].SetString(parameters["solver_settings"]["model_part_name"].GetString() + '.Results')
        from_json_check_result_settings["input_file_name"].SetString(results_filename)
        self.from_json_check_result = FromJsonCheckResultProcess(model, from_json_check_result_settings)
        self.from_json_check_result.ExecuteInitialize()
        self.from_json_check_result.ExecuteFinalizeSolutionStep()

    def __post_process(self, main_model_part, post_type = "gid"):
        if post_type == "gid":
            self.gid_output = GiDOutputProcess(
                main_model_part,
                main_model_part.Name,
                KratosMultiphysics.Parameters(r"""
                {
                    "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results"       : ["DISTANCE"],
                    "gauss_point_results" : []
                    }
                }"""))

            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()
            self.gid_output.ExecuteInitializeSolutionStep()
            self.gid_output.PrintOutput()
            self.gid_output.ExecuteFinalizeSolutionStep()
            self.gid_output.ExecuteFinalize()

        elif post_type == "vtk":
            vtk_output_parameters = KratosMultiphysics.Parameters(r"""
            {
                "model_part_name": "",
                "extrapolate_gauss_points": false,
                "nodal_solution_step_data_variables" : ["DISTANCE"],
                "gauss_point_variables": []
            }""")
            vtk_output_parameters["model_part_name"].SetString(main_model_part.Name)
            self.vtk_output_process = VtkOutputProcess(
                main_model_part.GetModel(),
                vtk_output_parameters)

            self.vtk_output_process.ExecuteInitialize()
            self.vtk_output_process.ExecuteBeforeSolutionLoop()
            self.vtk_output_process.ExecuteInitializeSolutionStep()
            self.vtk_output_process.PrintOutput()
            self.vtk_output_process.ExecuteFinalizeSolutionStep()
            self.vtk_output_process.ExecuteFinalize()

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
