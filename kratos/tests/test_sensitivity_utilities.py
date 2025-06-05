import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

def CalculateAnalyticalValue(nodal_neighbour_map, common_entity_map, node, derivative_dimension, value_dimension):
    analytical_value = KratosMultiphysics.Matrix((len(nodal_neighbour_map[node.Id]) + 1) * derivative_dimension, value_dimension, 0.0)
    for i_d in range(derivative_dimension):
        for v_d in range(value_dimension):
            for n_elem_id in common_entity_map[node.Id][node.Id]:
                analytical_value[i_d, v_d] += n_elem_id * node.Id * (i_d + 1) * (v_d + 1)
    for i_n, n_node_id in enumerate(nodal_neighbour_map[node.Id]):
        for i_d in range(derivative_dimension):
            for v_d in range(value_dimension):
                for common_entity in common_entity_map[node.Id][n_node_id]:
                    analytical_value[(i_n + 1) * derivative_dimension + i_d, v_d] += common_entity * n_node_id * (i_d + 1) * (v_d + 1)

    return analytical_value

def CalculateEntityNeighbourMaps(entities):
    nodal_neighbour_map = {}
    entity_neighbour_map = {}
    for entity in entities:
        for node_i in entity.GetGeometry():
            if (nodal_neighbour_map.get(node_i.Id) is None):
                nodal_neighbour_map[node_i.Id] = []

            if (entity_neighbour_map.get(node_i.Id) is None):
                entity_neighbour_map[node_i.Id] = []

            entity_neighbour_map[node_i.Id].append(entity.Id)

            for node_j in entity.GetGeometry():
                if (node_i.Id != node_j.Id and not (node_j.Id in nodal_neighbour_map[node_i.Id])):
                    nodal_neighbour_map[node_i.Id].append(node_j.Id)

    return nodal_neighbour_map, entity_neighbour_map

def CalculateCommonEntityMap(nodes, entity_neighbour_map, nodal_neighbour_map):
    common_entity_map = {}
    for node in nodes:
        if (entity_neighbour_map.get(node.Id) is not None and nodal_neighbour_map.get(node.Id) is not None):
            common_entity_map[node.Id] = {}
            common_entity_map[node.Id][node.Id] = entity_neighbour_map[node.Id]
            for neighbour_node_id in nodal_neighbour_map[node.Id]:
                current_entity_ids = []
                for base_neighbour_entity_id in entity_neighbour_map[node.Id]:
                    for derivative_neighbour_entity_id in entity_neighbour_map[neighbour_node_id]:
                        if (base_neighbour_entity_id == derivative_neighbour_entity_id):
                            current_entity_ids.append(base_neighbour_entity_id)
                            break

                common_entity_map[node.Id][neighbour_node_id] = current_entity_ids
    return common_entity_map


class TestSensitivityUtilitiesTwoDimSymmetricalSquare(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, self.model_part.Nodes)
        KratosMultiphysics.VariableUtils().ClearNonHistoricalData(self.model_part.Nodes)
        KratosMultiphysics.VariableUtils().ClearNonHistoricalData(self.model_part.Elements)
        KratosMultiphysics.VariableUtils().ClearNonHistoricalData(self.model_part.Conditions)

    def test_AssignElementDerivativesToNodes(self):
        nodal_neighbour_map = {61: [58, 59, 60], 60: [53, 54, 56, 58, 61], 29: [18, 22, 31, 38, 45], 28: [21, 25, 32, 37], 27: [20, 26, 31, 36], 26: [15, 17, 20, 23, 27, 33, 34, 36], 25: [15, 16, 21, 23, 28, 33, 35, 37], 24: [19, 21, 30, 32], 23: [15, 25, 26, 33], 22: [18, 20, 29, 31], 21: [12, 13, 16, 19, 24, 25, 28, 32], 20: [11, 14, 17, 18, 22, 26, 27, 31], 19: [9, 13, 21, 24, 30], 18: [8, 14, 20, 22, 29], 17: [11, 15, 20, 26], 16: [12, 15, 21, 25], 15: [10, 11, 12, 16, 17, 23, 25, 26], 14: [8, 11, 18, 20], 1: [2, 3, 4], 2: [1, 3, 4, 5], 3: [1, 2, 5, 6, 8], 4: [1, 2, 5, 7, 9], 5: [2, 3, 4, 6, 7, 10, 11, 12], 6: [3, 5, 8, 11], 7: [4, 5, 9, 12], 8: [3, 6, 11, 14, 18], 9: [4, 7, 12, 13, 19], 10: [5, 11, 12, 15], 11: [5, 6, 8, 10, 14, 15, 17, 20], 12: [5, 7, 9, 10, 13, 15, 16, 21], 13: [9, 12, 19, 21], 30: [19, 24, 32, 39, 46], 31: [20, 22, 27, 29, 36, 38, 40, 47], 32: [21, 24, 28, 30, 37, 39, 41, 48], 33: [23, 25, 26, 34, 35, 42, 43, 44], 34: [26, 33, 36, 43], 35: [25, 33, 37, 44], 36: [26, 27, 31, 34, 40, 43, 49, 52], 37: [25, 28, 32, 35, 41, 44, 50, 51], 38: [29, 31, 45, 47], 39: [30, 32, 46, 48], 40: [31, 36, 47, 52], 41: [32, 37, 48, 51], 42: [33, 43, 44, 53], 43: [33, 34, 36, 42, 49, 53, 54, 56], 44: [33, 35, 37, 42, 50, 53, 55, 57], 45: [29, 38, 47], 46: [30, 39, 48], 47: [31, 38, 40, 45, 52], 48: [32, 39, 41, 46, 51], 49: [36, 43, 52, 56], 50: [37, 44, 51, 57], 51: [37, 41, 48, 50, 57], 52: [36, 40, 47, 49, 56], 53: [42, 43, 44, 54, 55, 58, 59, 60], 54: [43, 53, 56, 60], 55: [44, 53, 57, 59], 56: [43, 49, 52, 54, 60], 57: [44, 50, 51, 55, 59], 58: [53, 59, 60, 61], 59: [53, 55, 57, 58, 61]}
        common_entity_map = {1: {1: [1, 4], 2: [1, 4], 3: [4], 4: [1]}, 2: {2: [1, 2, 3, 4], 1: [1, 4], 3: [3, 4], 4: [1, 2], 5: [2, 3]}, 3: {3: [3, 4, 21, 24], 1: [4], 2: [3, 4], 5: [3, 21], 6: [21, 24], 8: [24]}, 4: {4: [1, 2, 5, 8], 1: [1], 2: [1, 2], 5: [2, 8], 7: [5, 8], 9: [5]}, 5: {5: [2, 3, 7, 8, 21, 22, 25, 28], 2: [2, 3], 3: [3, 21], 4: [2, 8], 6: [21, 22], 7: [7, 8], 10: [25, 28], 11: [22, 28], 12: [7, 25]}, 6: {6: [21, 22, 23, 24], 3: [21, 24], 5: [21, 22], 8: [23, 24], 11: [22, 23]}, 7: {7: [5, 6, 7, 8], 4: [5, 8], 5: [7, 8], 9: [5, 6], 12: [6, 7]}, 8: {8: [23, 24, 41, 44], 3: [24], 6: [23, 24], 11: [23, 41], 14: [41, 44], 18: [44]}, 9: {9: [5, 6, 9, 12], 4: [5], 7: [5, 6], 12: [6, 12], 13: [9, 12], 19: [9]}, 10: {10: [25, 26, 27, 28], 5: [25, 28], 11: [27, 28], 12: [25, 26], 15: [26, 27]}, 11: {11: [22, 23, 27, 28, 41, 42, 45, 48], 5: [22, 28], 6: [22, 23], 8: [23, 41], 10: [27, 28], 14: [41, 42], 15: [27, 45], 17: [45, 48], 20: [42, 48]}, 12: {12: [6, 7, 11, 12, 25, 26, 29, 32], 5: [7, 25], 7: [6, 7], 9: [6, 12], 10: [25, 26], 13: [11, 12], 15: [26, 32], 16: [29, 32], 21: [11, 29]}, 13: {13: [9, 10, 11, 12], 9: [9, 12], 12: [11, 12], 19: [9, 10], 21: [10, 11]}, 14: {14: [41, 42, 43, 44], 8: [41, 44], 11: [41, 42], 18: [43, 44], 20: [42, 43]}, 15: {15: [26, 27, 31, 32, 45, 46, 49, 52], 10: [26, 27], 11: [27, 45], 12: [26, 32], 16: [31, 32], 17: [45, 46], 23: [49, 52], 25: [31, 49], 26: [46, 52]}, 16: {16: [29, 30, 31, 32], 12: [29, 32], 15: [31, 32], 21: [29, 30], 25: [30, 31]}, 17: {17: [45, 46, 47, 48], 11: [45, 48], 15: [45, 46], 20: [47, 48], 26: [46, 47]}, 18: {18: [43, 44, 61, 64], 8: [44], 14: [43, 44], 20: [43, 61], 22: [61, 64], 29: [64]}, 19: {19: [9, 10, 13, 16], 9: [9], 13: [9, 10], 21: [10, 16], 24: [13, 16], 30: [13]}, 20: {20: [42, 43, 47, 48, 61, 62, 65, 68], 11: [42, 48], 14: [42, 43], 17: [47, 48], 18: [43, 61], 22: [61, 62], 26: [47, 65], 27: [65, 68], 31: [62, 68]}, 21: {21: [10, 11, 15, 16, 29, 30, 33, 36], 12: [11, 29], 13: [10, 11], 16: [29, 30], 19: [10, 16], 24: [15, 16], 25: [30, 36], 28: [33, 36], 32: [15, 33]}, 22: {22: [61, 62, 63, 64], 18: [61, 64], 20: [61, 62], 29: [63, 64], 31: [62, 63]}, 23: {23: [49, 50, 51, 52], 15: [49, 52], 25: [49, 50], 26: [51, 52], 33: [50, 51]}, 24: {24: [13, 14, 15, 16], 19: [13, 16], 21: [15, 16], 30: [13, 14], 32: [14, 15]}, 25: {25: [30, 31, 35, 36, 49, 50, 53, 56], 15: [31, 49], 16: [30, 31], 21: [30, 36], 23: [49, 50], 28: [35, 36], 33: [50, 56], 35: [53, 56], 37: [35, 53]}, 26: {26: [46, 47, 51, 52, 65, 66, 69, 72], 15: [46, 52], 17: [46, 47], 20: [47, 65], 23: [51, 52], 27: [65, 66], 33: [51, 69], 34: [69, 72], 36: [66, 72]}, 27: {27: [65, 66, 67, 68], 20: [65, 68], 26: [65, 66], 31: [67, 68], 36: [66, 67]}, 28: {28: [33, 34, 35, 36], 21: [33, 36], 25: [35, 36], 32: [33, 34], 37: [34, 35]}, 29: {29: [63, 64, 81, 84], 18: [64], 22: [63, 64], 31: [63, 81], 38: [81, 84], 45: [84]}, 30: {30: [13, 14, 17, 20], 19: [13], 24: [13, 14], 32: [14, 20], 39: [17, 20], 46: [17]}, 31: {31: [62, 63, 67, 68, 81, 82, 85, 88], 20: [62, 68], 22: [62, 63], 27: [67, 68], 29: [63, 81], 36: [67, 85], 38: [81, 82], 40: [85, 88], 47: [82, 88]}, 32: {32: [14, 15, 19, 20, 33, 34, 37, 40], 21: [15, 33], 24: [14, 15], 28: [33, 34], 30: [14, 20], 37: [34, 40], 39: [19, 20], 41: [37, 40], 48: [19, 37]}, 33: {33: [50, 51, 55, 56, 69, 70, 73, 76], 23: [50, 51], 25: [50, 56], 26: [51, 69], 34: [69, 70], 35: [55, 56], 42: [73, 76], 43: [70, 76], 44: [55, 73]}, 34: {34: [69, 70, 71, 72], 26: [69, 72], 33: [69, 70], 36: [71, 72], 43: [70, 71]}, 35: {35: [53, 54, 55, 56], 25: [53, 56], 33: [55, 56], 37: [53, 54], 44: [54, 55]}, 36: {36: [66, 67, 71, 72, 85, 86, 89, 92], 26: [66, 72], 27: [66, 67], 31: [67, 85], 34: [71, 72], 40: [85, 86], 43: [71, 89], 49: [89, 92], 52: [86, 92]}, 37: {37: [34, 35, 39, 40, 53, 54, 57, 60], 25: [35, 53], 28: [34, 35], 32: [34, 40], 35: [53, 54], 41: [39, 40], 44: [54, 60], 50: [57, 60], 51: [39, 57]}, 38: {38: [81, 82, 83, 84], 29: [81, 84], 31: [81, 82], 45: [83, 84], 47: [82, 83]}, 39: {39: [17, 18, 19, 20], 30: [17, 20], 32: [19, 20], 46: [17, 18], 48: [18, 19]}, 40: {40: [85, 86, 87, 88], 31: [85, 88], 36: [85, 86], 47: [87, 88], 52: [86, 87]}, 41: {41: [37, 38, 39, 40], 32: [37, 40], 37: [39, 40], 48: [37, 38], 51: [38, 39]}, 42: {42: [73, 74, 75, 76], 33: [73, 76], 43: [75, 76], 44: [73, 74], 53: [74, 75]}, 43: {43: [70, 71, 75, 76, 89, 90, 93, 96], 33: [70, 76], 34: [70, 71], 36: [71, 89], 42: [75, 76], 49: [89, 90], 53: [75, 93], 54: [93, 96], 56: [90, 96]}, 44: {44: [54, 55, 59, 60, 73, 74, 77, 80], 33: [55, 73], 35: [54, 55], 37: [54, 60], 42: [73, 74], 50: [59, 60], 53: [74, 80], 55: [77, 80], 57: [59, 77]}, 45: {45: [83, 84], 29: [84], 38: [83, 84], 47: [83]}, 46: {46: [17, 18], 30: [17], 39: [17, 18], 48: [18]}, 47: {47: [82, 83, 87, 88], 31: [82, 88], 38: [82, 83], 40: [87, 88], 45: [83], 52: [87]}, 48: {48: [18, 19, 37, 38], 32: [19, 37], 39: [18, 19], 41: [37, 38], 46: [18], 51: [38]}, 49: {49: [89, 90, 91, 92], 36: [89, 92], 43: [89, 90], 52: [91, 92], 56: [90, 91]}, 50: {50: [57, 58, 59, 60], 37: [57, 60], 44: [59, 60], 51: [57, 58], 57: [58, 59]}, 51: {51: [38, 39, 57, 58], 37: [39, 57], 41: [38, 39], 48: [38], 50: [57, 58], 57: [58]}, 52: {52: [86, 87, 91, 92], 36: [86, 92], 40: [86, 87], 47: [87], 49: [91, 92], 56: [91]}, 53: {53: [74, 75, 79, 80, 93, 94, 97, 100], 42: [74, 75], 43: [75, 93], 44: [74, 80], 54: [93, 94], 55: [79, 80], 58: [97, 100], 59: [79, 97], 60: [94, 100]}, 54: {54: [93, 94, 95, 96], 43: [93, 96], 53: [93, 94], 56: [95, 96], 60: [94, 95]}, 55: {55: [77, 78, 79, 80], 44: [77, 80], 53: [79, 80], 57: [77, 78], 59: [78, 79]}, 56: {56: [90, 91, 95, 96], 43: [90, 96], 49: [90, 91], 52: [91], 54: [95, 96], 60: [95]}, 57: {57: [58, 59, 77, 78], 44: [59, 77], 50: [58, 59], 51: [58], 55: [77, 78], 59: [78]}, 58: {58: [97, 98, 99, 100], 53: [97, 100], 59: [97, 98], 60: [99, 100], 61: [98, 99]}, 59: {59: [78, 79, 97, 98], 53: [79, 97], 55: [78, 79], 57: [78], 58: [97, 98], 61: [98]}, 60: {60: [94, 95, 99, 100], 53: [94, 100], 54: [94, 95], 56: [95], 58: [99, 100], 61: [99]}, 61: {61: [98, 99], 58: [98, 99], 59: [98], 60: [99]}}

        # folowing lines are used merely in the case given mdpa changes then to generage new nodal_neighbour_map, and common_entity_map
        # nodal_neighbour_map, element_neighbour_map = CalculateEntityNeighbourMaps(self.model_part.Elements)
        # common_entity_map = CalculateCommonEntityMap(self.model_part.Nodes, element_neighbour_map, nodal_neighbour_map)

        self.__InitializeAssignEntityDerivativesToNodesValues(self.model_part.Elements, 3)

        KratosMultiphysics.SensitivityUtilities.AssignElementDerivativesToNodes(
            self.model_part,
            self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
            KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE,
            nodal_neighbour_map,
            1.0/3.0,
            KratosMultiphysics.SLIP,
            True)

        self.__CheckAssignEntityDerivativesToNodesValues(nodal_neighbour_map, common_entity_map)

    def test_AssignConditionDerivativesToNodes(self):
        nodal_neighbour_map = {60: [56, 61], 56: [60, 52], 52: [56, 47], 47: [52, 45], 45: [47, 29], 18: [8, 29], 8: [18, 3], 3: [8, 1], 46: [48, 30], 48: [46, 51], 29: [45, 18], 61: [60, 59], 9: [19, 4], 19: [9, 30], 30: [46, 19], 57: [59, 51], 59: [57, 61], 1: [4, 3], 4: [1, 9], 51: [48, 57]}
        common_entity_map = {1: {1: [14, 17], 4: [14], 3: [17]}, 3: {3: [6, 17], 8: [6], 1: [17]}, 4: {4: [14, 19], 1: [14], 9: [19]}, 8: {8: [5, 6], 18: [5], 3: [6]}, 9: {9: [10, 19], 19: [10], 4: [19]}, 18: {18: [5, 18], 8: [5], 29: [18]}, 19: {19: [10, 20], 9: [10], 30: [20]}, 29: {29: [8, 18], 45: [8], 18: [18]}, 30: {30: [11, 20], 46: [11], 19: [20]}, 45: {45: [4, 8], 47: [4], 29: [8]}, 46: {46: [7, 11], 48: [7], 30: [11]}, 47: {47: [3, 4], 52: [3], 45: [4]}, 48: {48: [7, 15], 46: [7], 51: [15]}, 51: {51: [15, 16], 48: [15], 57: [16]}, 52: {52: [2, 3], 56: [2], 47: [3]}, 56: {56: [1, 2], 60: [1], 52: [2]}, 57: {57: [12, 16], 59: [12], 51: [16]}, 59: {59: [12, 13], 57: [12], 61: [13]}, 60: {60: [1, 9], 56: [1], 61: [9]}, 61: {61: [9, 13], 60: [9], 59: [13]}}

        # folowing lines are used merely in the case given mdpa changes then to generage new nodal_neighbour_map, and common_entity_map
        # nodal_neighbour_map, condition_neighbour_map = CalculateEntityNeighbourMaps(self.model_part.Conditions)
        # common_entity_map = CalculateCommonEntityMap(self.model_part.Nodes, condition_neighbour_map, nodal_neighbour_map)

        self.__InitializeAssignEntityDerivativesToNodesValues(self.model_part.Conditions, 2)

        KratosMultiphysics.SensitivityUtilities.AssignConditionDerivativesToNodes(
            self.model_part,
            self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
            KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE,
            nodal_neighbour_map,
            1.0/3.0,
            KratosMultiphysics.SLIP,
            True)

        self.__CheckAssignEntityDerivativesToNodesValues(nodal_neighbour_map, common_entity_map)

    def __InitializeAssignEntityDerivativesToNodesValues(self, entities, nodes_per_entity):
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, entities)

        derivative_dimension = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        value_dimension = derivative_dimension

        for entity in entities:
            m = KratosMultiphysics.Matrix(nodes_per_entity * derivative_dimension, value_dimension)
            for i, node in enumerate(entity.GetGeometry()):
                node.Set(KratosMultiphysics.SLIP, True)
                for d_i in range(derivative_dimension):
                    for v_i in range(value_dimension):
                        m[i * derivative_dimension + d_i, v_i] = node.Id * (d_i+1) * (v_i + 1) * entity.Id * 3.0
            entity.SetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE, m)

        self.model_part.GetCommunicator().SynchronizeOrNodalFlags(KratosMultiphysics.SLIP)

    def __CheckAssignEntityDerivativesToNodesValues(self, nodal_neighbour_map, common_entity_map):
        derivative_dimension = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        value_dimension = derivative_dimension

        for node in self.model_part.Nodes:
            value = node.GetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE)
            if (node.Is(KratosMultiphysics.SLIP)):
                analytical_value = CalculateAnalyticalValue(nodal_neighbour_map, common_entity_map, node, derivative_dimension, value_dimension)
            else:
                analytical_value = KratosMultiphysics.Matrix(1, 1, 0.0)
            self.assertMatrixAlmostEqual(value, analytical_value, 12)

class TestSensitivityUtilitiesGetSensitivityVariableName(KratosUnittest.TestCase):
    def test_GetSensitivityVariableName(self):
        import KratosMultiphysics.StructuralMechanicsApplication
        self.assertEqual(KratosMultiphysics.SensitivityUtilities.GetSensitivityVariableName(KratosMultiphysics.POISSON_RATIO), "POISSON_RATIO_SENSITIVITY")
        self.assertEqual(KratosMultiphysics.SensitivityUtilities.GetSensitivityVariableName(KratosMultiphysics.SHAPE_SENSITIVITY), "SHAPE_SENSITIVITY")
        self.assertEqual(KratosMultiphysics.SensitivityUtilities.GetSensitivityVariableName(KratosMultiphysics.SHAPE_SENSITIVITY_X), "SHAPE_SENSITIVITY_X")

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

