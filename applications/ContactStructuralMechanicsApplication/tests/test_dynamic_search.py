from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestDynamicSearch(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _dynamic_search_tests(self, input_filename, num_nodes):
        self.main_model_part = KratosMultiphysics.ModelPart("Structure")
        
        ## Creation of the Kratos model (build sub_model_parts or submeshes)
        self.StructureModel = {"Structure": self.main_model_part}
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        
        self.main_model_part.CloneTimeStep(1.01)
        
        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NORMAL_CONTACT_STRESS, ContactStructuralMechanicsApplication.WEIGHTED_GAP, self.main_model_part)

        if (self.main_model_part.HasSubModelPart("Contact")):
            interface_model_part = self.main_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = self.main_model_part.CreateSubModelPart("Contact")
        
        self.contact_model_part = self.main_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2")
        
        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.SLAVE, False)
        del(node)
        model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        for node in model_part_slave.Nodes:
            node.Set(KratosMultiphysics.SLAVE, True)
            # DEBUG
            #node.X -= 9.81 / 32.0
            #node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, -9.81 / 32.0)
            node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_X, -9.81)
        del(node)
        model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")
        for node in model_part_master.Nodes:
            node.Set(KratosMultiphysics.MASTER, True)
        del(node)
        
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 1 
        self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 0.5 
        
        for prop in self.main_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = 3 
            
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = 3.0e-1
        
        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
            
        Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(self.main_model_part)

        interface_parameters = KratosMultiphysics.Parameters("""{"simplify_geometry": false}""")
        Preprocess.GenerateInterfacePart3D(self.main_model_part, self.contact_model_part, interface_parameters)
            
        # We copy the conditions to the ContactSubModelPart
        for cond in self.contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)    
        del(cond)
        for node in self.contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)    
        del(node)

        # We compute NODAL_H that can be used in the search and some values computation
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.contact_model_part)
        self.find_nodal_h.Execute()

        # We initialize the conditions    
        alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(self.contact_model_part) 
        alm_init_var.Execute()

        search_parameters = KratosMultiphysics.Parameters("""
        {
            "search_factor"               : 3.5,
            "allocation_size"             : 1000,
            "check_gap"                   : "MappingCheck",
            "type_search"                 : "InRadius"
        }
        """)
        if (num_nodes == 3):
            contact_search = ContactStructuralMechanicsApplication.TreeContactSearch3D3N(self.main_model_part, search_parameters)
        else:
            contact_search = ContactStructuralMechanicsApplication.TreeContactSearch3D4N(self.main_model_part, search_parameters)
        
        # We initialize the search utility
        contact_search.CreatePointListMortar()
        contact_search.InitializeMortarConditions()
        contact_search.UpdateMortarConditions()

        ## DEBUG
        #self.__post_process()
                
        import from_json_check_result_process

        check_parameters = KratosMultiphysics.Parameters("""
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
        
        check = from_json_check_result_process.FromJsonCheckResultProcess(self.StructureModel, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()
        
        #import json_output_process
        
        #out_parameters = KratosMultiphysics.Parameters("""
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

        #out = json_output_process.JsonOutputProcess(self.StructureModel, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()
                
    def test_dynamic_search_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/integration_tests/test_double_curvature_integration_triangle"
        
        self._dynamic_search_tests(input_filename, 3)
        
    def test_dynamic_search_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/integration_tests/test_double_curvature_integration_quadrilateral"
        
        self._dynamic_search_tests(input_filename, 4)
        
    def __post_process(self):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(self.main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "MultiFileFlag": "SingleFile"
                                                },        
                                                "nodal_results"       : ["NORMAL","DISPLACEMENT","VELOCITY","ACCELERATION"],
                                                "nodal_nonhistorical_results": ["DELTA_COORDINATES","AUXILIAR_COORDINATES","NORMAL_GAP"],
                                                "nodal_flags_results": ["ACTIVE","SLAVE"]
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
