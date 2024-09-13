import os
import copy

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.FSIApplication import fsi_analysis
from KratosMultiphysics.kratos_utilities import DeleteTimeFiles
from KratosMultiphysics.kratos_utilities import DeleteFilesEndingWith

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class EmbeddedFsiTest(KratosUnittest.TestCase):

    def setUp(self):
        self.print_output = False
        self.reference_values_check = True
        self.reference_values_output = False
        self.check_relative_tolerance = 1.0e-2
        self.check_absolute_tolerance = 1.0e-3
        self.work_folder = "embedded_fsi_test"
        self.settings = "ProjectParameters.json"

    def tearDown(self):
        DeleteTimeFiles(self.work_folder)
        if not self.print_output:
            DeleteFilesEndingWith(self.work_folder, ".post.bin")
            DeleteFilesEndingWith(self.work_folder, ".post.lst")

    def testEmbeddedVolumetricMVQN(self):
        self.fluid_element_settings = KratosMultiphysics.Parameters("""{
            "element_type" : "embedded_weakly_compressible_navier_stokes",
            "is_slip": true,
            "slip_length": 1.0e6,
            "penalty_coefficient": 1.0e3,
            "dynamic_tau": 1.0
        }""")
        self.structure_filename = "embedded_fsi_test_structure_volumetric"
        self.structure_materials_filename = "StructuralMaterialsVolumetric.json"
        self.convergence_accelerator = "MVQN"
        self.check_variables_list = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "LINE_LOAD_X", "LINE_LOAD_Y"]

        with WorkFolderScope(self.work_folder):
            model = KratosMultiphysics.Model()
            parameter_file = open(self.settings, 'r')
            self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
            self._CustomizeProjectParameters()
            fsi_analysis.FsiAnalysis(model, self.parameters).Run()

    def testEmbeddedVolumetricIBQNMVQN(self):
        self.fluid_element_settings = KratosMultiphysics.Parameters("""{
            "element_type" : "embedded_weakly_compressible_navier_stokes",
            "is_slip": true,
            "slip_length": 1.0e6,
            "penalty_coefficient": 1.0e3,
            "dynamic_tau": 1.0
        }""")
        self.structure_filename = "embedded_fsi_test_structure_volumetric"
        self.structure_materials_filename = "StructuralMaterialsVolumetric.json"
        self.convergence_accelerator = "IBQN_MVQN"
        self.check_variables_list = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "LINE_LOAD_X", "LINE_LOAD_Y"]

        with WorkFolderScope(self.work_folder):
            model = KratosMultiphysics.Model()
            parameter_file = open(self.settings, 'r')
            self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
            self._CustomizeProjectParameters()
            fsi_analysis.FsiAnalysis(model, self.parameters).Run()

    def testEmbeddedThinMVQN(self):
        self.fluid_element_settings = KratosMultiphysics.Parameters("""{
            "element_type" : "embedded_weakly_compressible_navier_stokes_discontinuous",
            "is_slip": true,
            "slip_length": 1.0e6,
            "penalty_coefficient": 1.0e3,
            "dynamic_tau": 1.0
        }""")
        self.structure_filename = "embedded_fsi_test_structure_thin"
        self.structure_materials_filename = "StructuralMaterialsThin.json"
        self.convergence_accelerator = "MVQN"
        self.check_variables_list = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "LINE_LOAD_X", "LINE_LOAD_Y"]

        with WorkFolderScope(self.work_folder):
            model = KratosMultiphysics.Model()
            parameter_file = open(self.settings, 'r')
            self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
            self._CustomizeProjectParameters()
            fsi_analysis.FsiAnalysis(model, self.parameters).Run()

    def testEmbeddedThinIBQNMVQN(self):
        self.fluid_element_settings = KratosMultiphysics.Parameters("""{
            "element_type" : "embedded_weakly_compressible_navier_stokes_discontinuous",
            "is_slip": true,
            "slip_length": 1.0e6,
            "penalty_coefficient": 1.0e3,
            "dynamic_tau": 1.0
        }""")
        self.structure_filename = "embedded_fsi_test_structure_thin"
        self.structure_materials_filename = "StructuralMaterialsThin.json"
        self.convergence_accelerator = "IBQN_MVQN"
        self.check_variables_list = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "LINE_LOAD_X", "LINE_LOAD_Y"]

        with WorkFolderScope(self.work_folder):
            model = KratosMultiphysics.Model()
            parameter_file = open(self.settings, 'r')
            self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
            self._CustomizeProjectParameters()
            fsi_analysis.FsiAnalysis(model, self.parameters).Run()

    def _CustomizeProjectParameters(self):
        self.parameters["solver_settings"]["fluid_solver_settings"]["formulation"] = self.fluid_element_settings
        self.parameters["solver_settings"]["structure_solver_settings"]["model_import_settings"]["input_filename"].SetString(self.structure_filename)
        self.parameters["solver_settings"]["structure_solver_settings"]["material_import_settings"]["materials_filename"].SetString(self.structure_materials_filename)
        self.parameters["solver_settings"]["coupling_settings"]["coupling_strategy_settings"]["solver_type"].SetString(self.convergence_accelerator)

        if self.print_output:
            self._AddOutput()

        if self.reference_values_output:
            self._AddReferenceValuesOutput()

        if self.reference_values_check:
            self._AddReferenceValuesCheck()

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "Parameters"    : {
                "model_part_name"        : "TO_BE_DEFINED",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"       : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"          : "time",
                        "output_control_type" : "step",
                        "output_interval"     : 1,
                        "body_output"         : true,
                        "node_output"         : false,
                        "skin_output"         : false,
                        "plane_output"        : [],
                        "nodal_results"       : [],
                        "gauss_point_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")

        fluid_output_settings = copy.deepcopy(gid_output_settings)
        fluid_output_settings["Parameters"]["model_part_name"].SetString("FluidModelPart")
        fluid_output_settings["Parameters"]["output_name"].SetString(f"{self.structure_filename}_{self.convergence_accelerator}_fluid")
        fluid_output_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"].SetStringArray(["VELOCITY","PRESSURE"])
        self.parameters["output_processes"]["gid_output"].Append(fluid_output_settings)

        structure_output_settings = copy.deepcopy(gid_output_settings)
        structure_output_settings["Parameters"]["model_part_name"].SetString("Structure")
        structure_output_settings["Parameters"]["output_name"].SetString(f"{self.structure_filename}_{self.convergence_accelerator}_structure")
        structure_output_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"].SetStringArray(["DISPLACEMENT"])
        self.parameters["output_processes"]["gid_output"].Append(structure_output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : [],
                "output_file_name" : "TO_BE_DEFINED",
                "model_part_name"  : "Structure",
                "time_frequency"   : 0.01
            }
        }""")
        output_file_name = f"{self.structure_filename}_{self.convergence_accelerator}_results.json"
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        json_output_settings["Parameters"]["output_variables"].SetStringArray(self.check_variables_list)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : [],
                "input_file_name"      : "TO_BE_DEFINED",
                "model_part_name"      : "Structure",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 0.01
            }
        }""")
        input_file_name = f"{self.structure_filename}_{self.convergence_accelerator}_results.json"
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["check_variables"].SetStringArray(self.check_variables_list)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_absolute_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

if __name__ == '__main__':
    KratosUnittest.main()

