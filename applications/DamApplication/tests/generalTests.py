import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

# Importing analysis
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

class KratosDamGeneralTests(KratosUnittest.TestCase):
    def setUp(self):
        self.print_results = False

    def test_joint_elastic_cohesive_2d_normal(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_elastic_cohesive_2d_normal/joint_elastic_cohesive_2d_normal_results.json"
            parameters_filename = "joint_elastic_cohesive_2d_normal/joint_elastic_cohesive_2d_normal_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_joint_elastic_cohesive_2d_shear(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_elastic_cohesive_2d_shear/joint_elastic_cohesive_2d_shear_results.json"
            parameters_filename = "joint_elastic_cohesive_2d_shear/joint_elastic_cohesive_2d_shear_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_joint_isotropic_damage_cohesive_2d_normal(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_isotropic_damage_cohesive_2d_normal/joint_isotropic_damage_cohesive_2d_normal_results.json"
            parameters_filename = "joint_isotropic_damage_cohesive_2d_normal/joint_isotropic_damage_cohesive_2d_normal_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_joint_isotropic_damage_cohesive_2d_shear(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_isotropic_damage_cohesive_2d_shear/joint_isotropic_damage_cohesive_2d_shear_results.json"
            parameters_filename = "joint_isotropic_damage_cohesive_2d_shear/joint_isotropic_damage_cohesive_2d_shear_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)


    def test_joint_elastic_cohesive_3d_normal(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_elastic_cohesive_3d_normal/joint_elastic_cohesive_3d_normal_results.json"
            parameters_filename = "joint_elastic_cohesive_3d_normal/joint_elastic_cohesive_3d_normal_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_joint_elastic_cohesive_3d_shear(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_elastic_cohesive_3d_shear/joint_elastic_cohesive_3d_shear_results.json"
            parameters_filename = "joint_elastic_cohesive_3d_shear/joint_elastic_cohesive_3d_shear_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_joint_isotropic_damage_cohesive_3d_normal(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_isotropic_damage_cohesive_3d_normal/joint_isotropic_damage_cohesive_3d_normal_results.json"
            parameters_filename = "joint_isotropic_damage_cohesive_3d_normal/joint_isotropic_damage_cohesive_3d_normal_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_joint_isotropic_damage_cohesive_3d_shear(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "joint_isotropic_damage_cohesive_3d_shear/joint_isotropic_damage_cohesive_3d_shear_results.json"
            parameters_filename = "joint_isotropic_damage_cohesive_3d_shear/joint_isotropic_damage_cohesive_3d_shear_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_mechanical(model, results_filename, parameters)
            self.__check_results_mechanical(model, results_filename, parameters)

    def test_construction(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "construction/construction_results.json"
            parameters_filename = "construction/construction_parameters.json"
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = DamAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results_thermo_mechanic(model, results_filename, parameters)
            self.__check_results_thermo_mechanic(model, results_filename, parameters)

    def __print_results_mechanical(self, model, results_filename, parameters):
        json_output_settings = KratosMultiphysics.Parameters(r"""
        {
            "output_variables": ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
            "output_file_name": "",
            "model_part_name": "",
            "time_frequency": 0.0
        }""")
        json_output_settings["model_part_name"].SetString(parameters["problem_data"]["model_part_name"].GetString() + '.Results')
        json_output_settings["output_file_name"].SetString(results_filename)
        self.json_output = JsonOutputProcess(model, json_output_settings)
        self.json_output.ExecuteInitialize()
        self.json_output.ExecuteBeforeSolutionLoop()
        self.json_output.ExecuteFinalizeSolutionStep()

    def __print_results_thermo_mechanic(self, model, results_filename, parameters):
        json_output_settings = KratosMultiphysics.Parameters(r"""
        {
            "output_variables": ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "TEMPERATURE"],
            "output_file_name": "",
            "model_part_name": "",
            "time_frequency": 0.0
        }""")
        json_output_settings["model_part_name"].SetString(parameters["problem_data"]["model_part_name"].GetString() + '.Results')
        json_output_settings["output_file_name"].SetString(results_filename)
        self.json_output = JsonOutputProcess(model, json_output_settings)
        self.json_output.ExecuteInitialize()
        self.json_output.ExecuteBeforeSolutionLoop()
        self.json_output.ExecuteFinalizeSolutionStep()

    def __check_results_mechanical(self, model, results_filename, parameters):
        from_json_check_result_settings = KratosMultiphysics.Parameters(r"""
        {
            "check_variables": ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
            "input_file_name": "",
            "model_part_name": ""
        }""")
        from_json_check_result_settings["model_part_name"].SetString(parameters["problem_data"]["model_part_name"].GetString() + '.Results')
        from_json_check_result_settings["input_file_name"].SetString(results_filename)
        self.from_json_check_result = FromJsonCheckResultProcess(model, from_json_check_result_settings)
        self.from_json_check_result.ExecuteInitialize()
        self.from_json_check_result.ExecuteFinalizeSolutionStep()

    def __check_results_thermo_mechanic(self, model, results_filename, parameters):
        from_json_check_result_settings = KratosMultiphysics.Parameters(r"""
        {
            "check_variables": ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "TEMPERATURE"],
            "input_file_name": "",
            "model_part_name": ""
        }""")
        from_json_check_result_settings["model_part_name"].SetString(parameters["problem_data"]["model_part_name"].GetString() + '.Results')
        from_json_check_result_settings["input_file_name"].SetString(results_filename)
        self.from_json_check_result = FromJsonCheckResultProcess(model, from_json_check_result_settings)
        self.from_json_check_result.ExecuteInitialize()
        self.from_json_check_result.ExecuteFinalizeSolutionStep()