from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Some imports
from KratosMultiphysics import from_json_check_result_process
#from KratosMultiphysics import json_output_process
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestCheckNormals(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _normal_check_process_tests(self, input_filename, custom_submodel_part = ""):
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        self.model = KM.Model()
        self.main_model_part = self.model.CreateModelPart("Main", 2)

        self.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)

        self.main_model_part.CloneTimeStep(1.01)

        KM.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        if custom_submodel_part == "":
            KM.VariableUtils().SetFlag(KM.INTERFACE, True, self.main_model_part.GetSubModelPart("CONTACT_Contact_slave_Auto1").Nodes)
            KM.VariableUtils().SetFlag(KM.INTERFACE, True, self.main_model_part.GetSubModelPart("CONTACT_Contact_master_Auto1").Nodes)
        else:
            KM.VariableUtils().SetFlag(KM.INTERFACE, True, self.main_model_part.GetSubModelPart(custom_submodel_part).Nodes)

        ## DEBUG
        #KM.ComputeNodesMeanNormalModelPart(self.main_model_part, True)

        # Check normals
        check_process = CSMA.NormalCheckProcess(self.main_model_part)
        check_process.Execute()

        ## DEBUG
        #self.__post_process()

        check_parameters = KM.Parameters("""
        {
            "check_variables"      : ["NORMAL"],
            "input_file_name"      : "",
            "model_part_name"      : "Main",
            "check_for_flag"       : "INTERFACE",
            "historical_value"     : true,
            "time_frequency"       : 0.0
        }
        """)

        check_parameters["input_file_name"].SetString(input_filename + "_check_normal.json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(self.model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        #out_parameters = KM.Parameters("""
        #{
            #"output_variables"     : ["NORMAL"],
            #"output_file_name"     : "",
            #"model_part_name"      : "Main",
            #"check_for_flag"       : "INTERFACE",
            #"historical_value"     : true,
            #"time_frequency"       : 0.0
        #}
        #""")

        #out_parameters["output_file_name"].SetString(input_filename + "_check_normal.json")

        #out = json_output_process.JsonOutputProcess(self.model, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_check_normals(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/inverted_normals"

        self._normal_check_process_tests(input_filename)

    def test_check_normals_quads(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/inverted_normals_quads"

        self._normal_check_process_tests(input_filename)

    def test_check_normals_s_shape(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/inverted_normals_s_shape"

        self._normal_check_process_tests(input_filename, "GENERIC_Contact_Auto1")

    def __post_process(self, debug = "GiD"):
        if debug == "GiD":
            self.gid_output = GiDOutputProcess(self.main_model_part,
                                        "gid_output",
                                        KM.Parameters("""
                                            {
                                                "result_file_configuration" : {
                                                    "gidpost_flags": {
                                                        "GiDPostMode": "GiD_PostBinary",
                                                        "WriteDeformedMeshFlag": "WriteUndeformed",
                                                        "MultiFileFlag": "SingleFile"
                                                    },
                                                    "nodal_results"       : ["NORMAL"],
                                                    "nodal_nonhistorical_results": [],
                                                    "nodal_flags_results": ["INTERFACE"]
                                                }
                                            }
                                            """)
                                        )

            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()
            self.gid_output.ExecuteInitializeSolutionStep()
            self.gid_output.PrintOutput()
            self.gid_output.ExecuteFinalizeSolutionStep()
            self.gid_output.ExecuteFinalize()
        elif debug == "VTK":
            self.vtk_output_process = VtkOutputProcess(self.model,
                                        KM.Parameters("""{
                                                "model_part_name"                    : "Main",
                                                "nodal_solution_step_data_variables" : ["NORMAL"],
                                                "nodal_data_value_variables": [],
                                                "nodal_flags" : ["INTERFACE"]
                                            }
                                            """)
                                        )

            self.vtk_output_process.ExecuteInitialize()
            self.vtk_output_process.ExecuteBeforeSolutionLoop()
            self.vtk_output_process.ExecuteInitializeSolutionStep()
            self.vtk_output_process.PrintOutput()
            self.vtk_output_process.ExecuteFinalizeSolutionStep()
            self.vtk_output_process.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
