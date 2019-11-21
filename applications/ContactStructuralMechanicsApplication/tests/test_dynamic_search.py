from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Some imports
from KratosMultiphysics import from_json_check_result_process
#from KratosMultiphysics import json_output_process
from KratosMultiphysics.gid_output_process import GiDOutputProcess

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestDynamicSearch(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _dynamic_search_tests(self, input_filename, num_nodes):
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        self.model = KM.Model()
        self.main_model_part = self.model.CreateModelPart("Structure", 2)

        self.main_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_H)

        self.main_model_part.CloneTimeStep(1.01)

        KM.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X,self.main_model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y,self.main_model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z,self.main_model_part)
        KM.VariableUtils().AddDof(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, CSMA.WEIGHTED_GAP, self.main_model_part)

        if self.main_model_part.HasSubModelPart("Contact"):
            interface_model_part = self.main_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = self.main_model_part.CreateSubModelPart("Contact")

        self.contact_model_part = self.main_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2")

        model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        #model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")
        KM.VariableUtils().SetFlag(KM.SLAVE, False, self.contact_model_part.Nodes)
        KM.VariableUtils().SetFlag(KM.MASTER, True, self.contact_model_part.Nodes)
        KM.VariableUtils().SetFlag(KM.SLAVE, True, model_part_slave.Nodes)
        KM.VariableUtils().SetFlag(KM.MASTER, False, model_part_slave.Nodes)

        for node in model_part_slave.Nodes:
            # DEBUG
            #node.X -= 9.81 / 32.0
            #node.SetSolutionStepValue(KM.DISPLACEMENT_X, -9.81 / 32.0)
            node.SetSolutionStepValue(KM.ACCELERATION_X, 1, -9.81)

        self.main_model_part.ProcessInfo[KM.STEP] = 1
        self.main_model_part.ProcessInfo[KM.DELTA_TIME] = 0.5

        for prop in self.main_model_part.GetProperties():
            prop[CSMA.INTEGRATION_ORDER_CONTACT] = 3

        self.main_model_part.ProcessInfo[CSMA.ACTIVE_CHECK_FACTOR] = 3.0e-1

        KM.VariableUtils().SetFlag(KM.INTERFACE, True, self.contact_model_part.Nodes)

        pre_process = CSMA.InterfacePreprocessCondition(self.main_model_part)

        interface_parameters = KM.Parameters("""{"simplify_geometry": false}""")
        pre_process.GenerateInterfacePart(self.contact_model_part, interface_parameters)

        # We copy the conditions to the ContactSubModelPart
        for cond in self.contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)
        for node in self.contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)

        # We compute NODAL_H that can be used in the search and some values computation
        self.find_nodal_h = KM.FindNodalHProcess(self.contact_model_part)
        self.find_nodal_h.Execute()

        # We initialize the conditions
        alm_init_var = CSMA.ALMFastInit(self.contact_model_part)
        alm_init_var.Execute()

        search_parameters = KM.Parameters("""
        {
            "dynamic_search"               : true,
            "simple_search"                : false,
            "normal_orientation_threshold" : 0.0
        }
        """)
        contact_search = CSMA.ContactSearchProcess(self.main_model_part, search_parameters)

        # We initialize the search utility
        contact_search.ExecuteInitialize()
        contact_search.ExecuteInitializeSolutionStep()

        ## DEBUG
        #self.__post_process()

        check_parameters = KM.Parameters("""
        {
            "check_variables"      : ["NORMAL_GAP"],
            "input_file_name"      : "",
            "model_part_name"      : "Structure",
            "historical_value"     : false,
            "time_frequency"       : 0.0,
            "sub_model_part_name"  : "Parts_Parts_Auto1"
        }
        """)

        check_parameters["input_file_name"].SetString(input_filename + "_dynamic_search.json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(self.model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        #out_parameters = KM.Parameters("""
        #{
            #"output_variables"     : ["NORMAL_GAP"],
            #"output_file_name"     : "",
            #"model_part_name"      : "Structure",
            #"historical_value"     : false,
            #"time_frequency"       : 0.0,
            #"sub_model_part_name"  : "Parts_Parts_Auto1"
        #}
        #""")

        #out_parameters["output_file_name"].SetString(input_filename + "_dynamic_search.json")

        #out = json_output_process.JsonOutputProcess(self.model, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_dynamic_search_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/integration_tests/test_double_curvature_integration_triangle"

        self._dynamic_search_tests(input_filename, 3)

    def test_dynamic_search_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/integration_tests/test_double_curvature_integration_quadrilateral"

        self._dynamic_search_tests(input_filename, 4)

    def __post_process(self):
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
                                                "nodal_results"       : ["NORMAL","DISPLACEMENT","VELOCITY","ACCELERATION"],
                                                "nodal_nonhistorical_results": ["DELTA_COORDINATES","AUXILIAR_COORDINATES","NORMAL_GAP"],
                                                "nodal_flags_results": ["ACTIVE","SLAVE","MASTER"]
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

if __name__ == '__main__':
    KratosUnittest.main()
