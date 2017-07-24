from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestDoubleCurvatureIntegration(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
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

        # We initialize the conditions    
        alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(contact_model_part) 
        alm_init_var.Execute()

        search_parameters = KratosMultiphysics.Parameters("""
        {
            "search_factor"               : 2.5,
            "allocation_size"             : 1000,
            "type_search"                 : "InRadius",
            "use_exact_integration"       : true
        }
        """)
        contact_search = ContactStructuralMechanicsApplication.TreeContactSearch(main_model_part, search_parameters)
        
        # We initialize the search utility
        contact_search.CreatePointListMortar()
        contact_search.InitializeMortarConditions()
        contact_search.UpdateMortarConditions()

        ## DEBUG
        #self.__post_process(main_model_part)
        #exact_integration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N(3, True)
        
        exact_integration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N(3)
        
        # These conditions are in the border, and can not be integrated 100% accurate
        list_of_border_cond = [1262,1263,1264,1265,1269,1270,1273,1275,1278,1282,1284,1285,1286,1288,1290,1291,1292,1294,1295,1297,1298,1302,1303,1305,1306,1307,1310,1313,1314,1318,1319,1320,1323,1325,1327,1328,1329,1331,1336,1337,1338,1340,1341,1342,1343,1344,1346,1347,1348,1349,1350,1353,1355,1357,1359,1360,1366,1367,1368,1369,1370,1377,1378,1379,1381,1382,1384,1385,1387,1393,1394,1395,1399,1400,1406,1410,1411,1412,1414,1415,1418,1419,1420,1424,1427,1429,1431,1436,1438,1444,1446,1447,1448,1449,1459,1462,1463,1465,1467,1468,1474,1477,1479,1485,1491,1493,1507,1515,1517,1531,1537,1539,1547,1549,1553,1563,1569,1575,1623,1640,1644,1654,1656,1663,1667,1675,1685,1687,1693,1697,1703,1707,1713,1715,1717,1719,1721,1723,1725]
        
        for cond in contact_model_part.Conditions:
            if cond.Is(KratosMultiphysics.SLAVE):
                if (cond.Id in list_of_border_cond == False):
                    area = exact_integration.TestGetExactAreaIntegration(cond)
                    condition_area = cond.GetArea()
                    self.assertAlmostEqual(area, condition_area)
                
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
                
                
                
                

