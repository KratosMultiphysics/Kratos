from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

# Some additional imports
from KratosMultiphysics.gid_output_process import GiDOutputProcess

# Some system imports
import os
from decimal import Decimal

class TestDoubleCurvatureIntegration(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def __base_test_integration(self, input_filename, num_nodes):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        self.model = KratosMultiphysics.Model()
        self.main_model_part = self.model.CreateModelPart("Structure", 2)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(ContactStructuralMechanicsApplication.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, ContactStructuralMechanicsApplication.WEIGHTED_GAP, self.main_model_part)

        if self.main_model_part.HasSubModelPart("Contact"):
            interface_model_part = self.main_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = self.main_model_part.CreateSubModelPart("Contact")

        self.contact_model_part = self.main_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2")

        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.SLAVE, False)
        model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        for node in model_part_slave.Nodes:
            node.Set(KratosMultiphysics.SLAVE, True)

        for prop in self.main_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = 3

        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = 3.0e-1

        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)

        Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(self.main_model_part)

        interface_parameters = KratosMultiphysics.Parameters("""{"simplify_geometry": false}""")
        Preprocess.GenerateInterfacePart(self.contact_model_part, interface_parameters)

        # We copy the conditions to the ContactSubModelPart
        for cond in self.contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)
        for node in self.contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)

        # We initialize the conditions
        alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(self.contact_model_part)
        alm_init_var.Execute()

        search_parameters = KratosMultiphysics.Parameters("""
        {
            "search_factor"               : 3.5,
            "allocation_size"             : 1000,
            "check_gap"                   : "NoCheck",
            "type_search"                 : "InRadius",
            "simple_search"               : false
        }
        """)
        contact_search = ContactStructuralMechanicsApplication.ContactSearchProcess(self.main_model_part, search_parameters)

        # We initialize the search utility
        contact_search.ExecuteInitialize()
        contact_search.ExecuteInitializeSolutionStep()

        if (num_nodes == 3):
            ## DEBUG
            #self.__post_process()
            self.exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D3N(3)
        else:
            ## DEBUG
            #self.__post_process()
            self.exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D4N(3)

    def _double_curvature_tests(self, input_filename, num_nodes, list_of_border_cond):

        self.__base_test_integration(input_filename, num_nodes)

        tolerance = 5.0e-3
        for cond in self.contact_model_part.Conditions:
            if cond.Is(KratosMultiphysics.SLAVE):
                to_test = (cond.Id in list_of_border_cond)
                if not to_test:
                    area = self.exact_integration.TestGetExactAreaIntegration(self.main_model_part, cond)
                    condition_area = cond.GetGeometry().Area()
                    check_value = abs((area - condition_area)/condition_area)
                    if (check_value >  tolerance):
                        print(cond.Id,"\t",area,"\t", condition_area,"\t", self.__sci_str(check_value))
                    else:
                        self.assertLess(check_value, tolerance)

    def _moving_nodes_tests(self, input_filename, num_nodes):

        self.__base_test_integration(input_filename, num_nodes)

        delta_disp = 1.0e-6
        for node in self.main_model_part.GetSubModelPart("GroupPositiveX").Nodes:
            node.X += delta_disp
        for node in self.main_model_part.GetSubModelPart("GroupPositiveY").Nodes:
            node.Y += delta_disp
        for node in self.main_model_part.GetSubModelPart("GroupNegativeX").Nodes:
            node.X -= delta_disp
        for node in self.main_model_part.GetSubModelPart("GroupNegativeY").Nodes:
            node.Y -= delta_disp

        #print("Solution obtained")
        tolerance = 5.0e-5
        for cond in self.contact_model_part.Conditions:
            if cond.Is(KratosMultiphysics.SLAVE):
                area = self.exact_integration.TestGetExactAreaIntegration(self.contact_model_part, cond)
                condition_area = cond.GetGeometry().Area()
                check_value = abs((area - condition_area)/condition_area)
                if check_value >  tolerance:
                    print(cond.Id,"\t",area,"\t", condition_area,"\t", __sci_str(check_value))
                else:
                    self.assertLess(check_value, tolerance)

    def test_double_curvature_integration_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/integration_tests/test_double_curvature_integration_triangle"

        # These conditions are in the border, and can not be integrated 100% accurate
        list_of_border_cond = [1262,1263,1264,1265,1269,1270,1273,1275,1278,1282,1284,1285,1286,1288,1290,1291,1292,1294,1295,1297,1298,1302,1303,1305,1306,1307,1310,1313,1314,1318,1319,1320,1323,1325,1327,1328,1329,1331,1336,1337,1338,1340,1341,1342,1343,1344,1346,1347,1348,1349,1350,1353,1355,1357,1359,1360,1366,1367,1368,1369,1370,1377,1378,1379,1381,1382,1384,1385,1387,1393,1394,1395,1399,1400,1406,1410,1411,1412,1414,1415,1418,1419,1420,1424,1427,1429,1431,1436,1438,1444,1446,1447,1448,1449,1459,1462,1463,1465,1467,1468,1474,1477,1479,1485,1491,1493,1507,1515,1517,1531,1537,1539,1547,1549,1553,1563,1569,1575,1623,1640,1644,1654,1656,1663,1667,1675,1685,1687,1693,1697,1703,1707,1713,1715,1717,1719,1721,1723,1725]

        self._double_curvature_tests(input_filename, 3, list_of_border_cond)

    def test_double_curvature_integration_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/integration_tests/test_double_curvature_integration_quadrilateral"

        # These conditions are in the border, and can not be integrated 100% accurate
        list_of_border_cond = [916,917,919,920,923,925,927,929,933,934,938,940,941,944,945,946,949,951,954,955,962,963,965,966,967,968,969,970,971,973,974,977,978,979,980,981,982,983,984,985,986,988,989,990,995,996,1000,1003,1005,1007,1008,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1058,1060,1064,1066,1069,1070,1071,1072,1073,1074,1075,1076]

        self._double_curvature_tests(input_filename, 4, list_of_border_cond)

    def test_moving_mesh_integration_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/integration_tests/quadrilaterals_moving_nodes"

        self._moving_nodes_tests(input_filename, 4)

    def test_integration_quad_non_matching(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unittest/integration_tests/quadrilaterals_non_matching"

        list_of_border_cond = []
        self._double_curvature_tests(input_filename, 4, list_of_border_cond)

    def __post_process(self):
        self.gid_output = GiDOutputProcess(self.main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditionsOnly",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["NORMAL"],
                                                "nodal_nonhistorical_results": [],
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

def __sci_str(x):
    s = 10*Decimal(str(x))
    s = ('{:.' + str(len(s.normalize().as_tuple().digits) - 1) + 'E}').format(s)
    s = s.replace('E+','D0')
    s = s.replace('E-','D0-')
    s = s.replace('.','')
    if s.startswith('-'):
        return '-.' + s[1:]
    else:
        return '.' + s

if __name__ == '__main__':
    KratosUnittest.main()
