from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestDoubleCurvatureIntegration(KratosUnittest.TestCase):
    def test_double_curvature_integration(self):
        input_filename = "integration_tests/test_double_curvature_integration"

        main_model_part = KratosMultiphysics.ModelPart("Structure")
        
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)
        main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        
        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NORMAL_CONTACT_STRESS, ContactStructuralMechanicsApplication.WEIGHTED_GAP, main_model_part)

        if (main_model_part.HasSubModelPart("Contact")):
            interface_model_part = main_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = main_model_part.CreateSubModelPart("Contact")
        
        contact_model_part = main_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2")
        
        for node in contact_model_part.Nodes:
            node.Set(KratosMultiphysics.SLAVE, False)
        del(node)
        model_part_slave = main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        for node in model_part_slave.Nodes:
            node.Set(KratosMultiphysics.SLAVE, True)
        del(node)
        
        for prop in main_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = 3 
            prop[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = 3.0e-1
        
        for node in contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
            
        Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(main_model_part)

        interface_parameters = KratosMultiphysics.Parameters("""{"condition_name": "", "final_string": "", "simplify_geometry": false}""")
        interface_parameters["condition_name"].SetString("ALMFrictionlessMortarContact")
        Preprocess.GenerateInterfacePart3D(main_model_part, contact_model_part, interface_parameters)
            
        # We copy the conditions to the ContactSubModelPart
        for cond in contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)    
        del(cond)
        for node in contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)    
        del(node)

        search_parameters = KratosMultiphysics.Parameters("""
        {
            "search_factor"               : 2.5,
            "allocation_size"             : 1000,
            "type_search"                 : "InRadius",
            "use_exact_integration"       : true
        }
        """)
        contact_search = ContactStructuralMechanicsApplication.TreeContactSearch(main_model_part, search_parameters)

        # We initialize the conditions    
        alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(contact_model_part) 
        alm_init_var.Execute()
        
        # We initialize the search utility
        contact_search.CreatePointListMortar()
        contact_search.InitializeMortarConditions()
        contact_search.UpdateMortarConditions()

        exact_integration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N(3)
        
        for cond in contact_model_part.Conditions:
            if cond.Is(KratosMultiphysics.SLAVE):
                area = 0
                exact_integration.TestGetExactAreaIntegration(cond, area)
                self.assertAlmostEqual(area, cond.GetArea)
                
    def __post_process(self, main_model_part):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },        
                                                "nodal_results"       : ["DISPLACEMENT","NORMAL_CONTACT_STRESS","WEIGHTED_GAP"],
                                                "nodal_nonhistorical_results": ["NORMAL","AUGMENTED_NORMAL_CONTACT_PRESSURE"],
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
                
                
                
                

