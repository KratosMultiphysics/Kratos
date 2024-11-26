import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing analysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Importing post-process
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

class TestCookMembrane(KratosUnittest.TestCase):
    def setUp(self):
        self.print_output = False
        self.print_results = False

    def test_cook_membrane_compressible_asgs(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "cook_membrane_test/cook_membrane_compressible_asgs"
            parameters_filename = "cook_membrane_test/cook_membrane_parameters.json"
            check_variables = ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"]

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Set materials settings
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("cook_membrane_test/cook_membrane_compressible_materials.json")

            # Set ASGS element settings
            parameters["modelers"][1]["parameters"]["elements_list"][0]["element_name"].SetString("SmallDisplacementMixedVolumetricStrainElement2D3N")

            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results(model, f"{results_filename}_results.json", check_variables)
            if self.print_output:
                output_variables = ["DISPLACEMENT", "VOLUMETRIC_STRAIN"]
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()), results_filename, output_variables)
            self.__check_results(model, f"{results_filename}_results.json", check_variables)

    def test_cook_membrane_incompressible_asgs(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "cook_membrane_test/cook_membrane_incompressible_asgs"
            parameters_filename = "cook_membrane_test/cook_membrane_parameters.json"
            check_variables = ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"]

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Set materials settings
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("cook_membrane_test/cook_membrane_incompressible_materials.json")

            # Set ASGS element settings
            parameters["modelers"][1]["parameters"]["elements_list"][0]["element_name"].SetString("SmallDisplacementMixedVolumetricStrainElement2D3N")

            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results(model, f"{results_filename}_results.json", check_variables)
            if self.print_output:
                output_variables = ["DISPLACEMENT", "VOLUMETRIC_STRAIN"]
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()), results_filename, output_variables)
            self.__check_results(model, f"{results_filename}_results.json", check_variables)

    def test_cook_membrane_incompressible_oss(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "cook_membrane_test/cook_membrane_incompressible_oss"
            parameters_filename = "cook_membrane_test/cook_membrane_parameters.json"
            check_variables = ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION_X","DISPLACEMENT_PROJECTION_Y", "VOLUMETRIC_STRAIN_PROJECTION"]

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Set materials settings
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("cook_membrane_test/cook_membrane_incompressible_materials.json")

            # Set linearised OSS element settings
            parameters["solver_settings"]["analysis_type"].SetString("non_linear") # Note that the problem becomes inherently non linear as we are linearising the projections
            parameters["solver_settings"]["use_orthogonal_subscales"].SetBool(True) # Note that this flag activates the residual projections for the linearised OSS
            parameters["modelers"][1]["parameters"]["elements_list"][0]["element_name"].SetString("SmallDisplacementMixedVolumetricStrainOssElement2D3N")

            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results(model, f"{results_filename}_results.json", check_variables)
            if self.print_output:
                output_variables = ["DISPLACEMENT", "VOLUMETRIC_STRAIN", "DISPLACEMENT_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION"]
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()), results_filename, output_variables)
            self.__check_results(model, f"{results_filename}_results.json", check_variables)

    def test_cook_membrane_incompressible_oss_non_linear(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "cook_membrane_test/cook_membrane_incompressible_oss_non_linear"
            parameters_filename = "cook_membrane_test/cook_membrane_parameters.json"
            check_variables = ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION_X","DISPLACEMENT_PROJECTION_Y", "VOLUMETRIC_STRAIN_PROJECTION"]

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Set materials settings
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("cook_membrane_test/cook_membrane_incompressible_materials.json")

            # Set non linear OSS element settings
            parameters["solver_settings"]["analysis_type"].SetString("non_linear") # Note that the problem is necessarily non linear
            parameters["solver_settings"]["use_orthogonal_subscales"].SetBool(False) # Note that the flag is not needed in here as the projections are actually solved in the linear system
            parameters["modelers"][1]["parameters"]["elements_list"][0]["element_name"].SetString("SmallDisplacementMixedVolumetricStrainOssNonLinearElement2D3N")

            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results(model, f"{results_filename}_results.json", check_variables)
            if self.print_output:
                output_variables = ["DISPLACEMENT", "VOLUMETRIC_STRAIN", "DISPLACEMENT_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION"]
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()), results_filename, output_variables)
            self.__check_results(model, f"{results_filename}_results.json", check_variables)

    def test_cook_membrane_incompressible_asgs_dynamic(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "cook_membrane_test/cook_membrane_incompressible_asgs_dynamic"
            parameters_filename = "cook_membrane_test/cook_membrane_parameters.json"
            check_variables = ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"]

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Set materials settings
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("cook_membrane_test/cook_membrane_incompressible_materials.json")

            # Set ASGS element settings
            parameters["modelers"][1]["parameters"]["elements_list"][0]["element_name"].SetString("SmallDisplacementMixedVolumetricStrainElement2D3N")

            # Set dynamic settings
            parameters["problem_data"]["end_time"].SetDouble(3.0)
            parameters["solver_settings"]["solver_type"].SetString("dynamic")

            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results(model, f"{results_filename}_results.json", check_variables)
            if self.print_output:
                output_variables = ["DISPLACEMENT", "VOLUMETRIC_STRAIN"]
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()), results_filename, output_variables)
            self.__check_results(model, f"{results_filename}_results.json", check_variables)

    def test_cook_membrane_incompressible_oss_dynamic(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            results_filename = "cook_membrane_test/cook_membrane_incompressible_oss_dynamic"
            parameters_filename = "cook_membrane_test/cook_membrane_parameters.json"
            check_variables = ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION_X","DISPLACEMENT_PROJECTION_Y", "VOLUMETRIC_STRAIN_PROJECTION"]

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Set materials settings
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("cook_membrane_test/cook_membrane_incompressible_materials.json")

            # Set linearised OSS element settings
            parameters["solver_settings"]["analysis_type"].SetString("non_linear") # Note that the problem becomes inherently non linear as we are linearising the projections
            parameters["solver_settings"]["use_orthogonal_subscales"].SetBool(True) # Note that this flag activates the residual projections for the linearised OSS
            parameters["modelers"][1]["parameters"]["elements_list"][0]["element_name"].SetString("SmallDisplacementMixedVolumetricStrainOssElement2D3N")

            # Set dynamic settings
            parameters["problem_data"]["end_time"].SetDouble(3.0)
            parameters["solver_settings"]["solver_type"].SetString("dynamic")

            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsAnalysis(model, parameters)
            simulation.Run()

            if self.print_results:
                self.__print_results(model, f"{results_filename}_results.json", check_variables)
            if self.print_output:
                output_variables = ["DISPLACEMENT", "VOLUMETRIC_STRAIN", "DISPLACEMENT_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION"]
                self.__post_process(model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString()), results_filename, output_variables)
            self.__check_results(model, f"{results_filename}_results.json", check_variables)

    def __print_results(self, model, results_filename, output_variables):
        json_output_settings = KratosMultiphysics.Parameters(r"""
        {
            "output_variables": ["DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"],
            "output_file_name": "",
            "time_frequency": 0.00,
            "model_part_name": "cook_membrane.Parts_ResultsCheck"
        }""")
        json_output_settings["output_file_name"].SetString(results_filename)
        json_output_settings["output_variables"].SetStringArray(output_variables)
        self.json_output = JsonOutputProcess(
            model,
            json_output_settings)
        self.json_output.ExecuteInitialize()
        self.json_output.ExecuteBeforeSolutionLoop()
        self.json_output.ExecuteFinalizeSolutionStep()

    def __check_results(self, model, results_filename, check_variables):
        from_json_check_result_settings = KratosMultiphysics.Parameters(r"""
        {
            "check_variables": [],
            "input_file_name": "",
            "model_part_name": "cook_membrane.Parts_ResultsCheck"
        }""")
        from_json_check_result_settings["input_file_name"].SetString(results_filename)
        from_json_check_result_settings["check_variables"].SetStringArray(check_variables)
        self.from_json_check_result = FromJsonCheckResultProcess(
            model,
            from_json_check_result_settings)
        self.from_json_check_result.ExecuteInitialize()
        self.from_json_check_result.ExecuteFinalizeSolutionStep()

    def __post_process(self, main_model_part, results_filename, nodal_results_variables):
        gid_output_parameters = KratosMultiphysics.Parameters(r"""{
            "result_file_configuration" : {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostBinary",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "nodal_results"       : [],
                "gauss_point_results" : []
            }
        }""")
        gid_output_parameters["result_file_configuration"]["nodal_results"].SetStringArray(nodal_results_variables)
        self.gid_output = GiDOutputProcess(
            main_model_part,
            results_filename,
            gid_output_parameters)

        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
