import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSModellingApplication as KratosRANS
from process_factory import KratosProcessFactory

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import random, math


class CustomProcessTest(UnitTest.TestCase):
    def testApplyFlagProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.submodelpart_1",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "STRUCTURE",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : "all"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.submodelpart_2",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "INLET",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : "all"
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        check_values = [True, True, True, False, False, False]
        for node, value in zip(self.model_part.Nodes, check_values):
            self.assertEqual(node.Is(Kratos.STRUCTURE), value)

        check_values = [False, False, True, True, True, False]
        for node, value in zip(self.model_part.Nodes, check_values):
            self.assertEqual(node.Is(Kratos.INLET), value)

        check_values = [True, True, False, False, False, False]
        for condition, value in zip(self.model_part.Conditions, check_values):
            self.assertEqual(condition.Is(Kratos.STRUCTURE), value)

        check_values = [False, False, True, True, False, False]
        for condition, value in zip(self.model_part.Conditions, check_values):
            self.assertEqual(condition.Is(Kratos.INLET), value)

    def testScalarCellCenterAveragingProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "FindNodalNeighboursProcess",
                "Parameters" : {
                    "model_part_name"      : "test"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "FindConditionParentProcess",
                "Parameters" : {
                    "model_part_name"      : "test.submodelpart_1"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyScalarCellCenteredAveragingProcess",
                "Parameters" : {
                    "echo_level"           : 0,
                    "model_part_name"      : "test.submodelpart_1",
                    "input_variable_name"  : "DENSITY",
                    "output_variable_name" : "DENSITY"
                }
            }
        ]''')

        check_values = []
        for node in self.model_part.Nodes:
            check_values.append(node.GetSolutionStepValue(Kratos.DENSITY))

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        check_values[0] = check_values[5] / 3.0
        check_values[1] = check_values[5] / 6.0 + check_values[4] / 6.0
        check_values[2] = check_values[4] / 3.0

        for node, value in zip(self.model_part.Nodes, check_values):
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.DENSITY), value, 9)

    def testVectorCellCenterAveragingProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "FindNodalNeighboursProcess",
                "Parameters" : {
                    "model_part_name"      : "test"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "FindConditionParentProcess",
                "Parameters" : {
                    "model_part_name"      : "test.submodelpart_1"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyVectorCellCenteredAveragingProcess",
                "Parameters" : {
                    "echo_level"           : 0,
                    "model_part_name"      : "test.submodelpart_1",
                    "input_variable_name"  : "VELOCITY",
                    "output_variable_name" : "VELOCITY"
                }
            }
        ]''')

        check_values = []
        for node in self.model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            check_values.append(
                [float(velocity[0]),
                 float(velocity[1]),
                 float(velocity[2])])

        for i in range(3):
            check_values[0][i] = check_values[5][i] / 3.0
            check_values[1][
                i] = check_values[5][i] / 6.0 + check_values[4][i] / 6.0
            check_values[2][i] = check_values[4][i] / 3.0

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node, value in zip(self.model_part.Nodes, check_values):
            node_value = node.GetSolutionStepValue(Kratos.VELOCITY)
            for i in range(3):
                self.assertAlmostEqual(node_value[i], value[i], 9)

    def testVectorAlignProcessTangential(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyVectorAlignProcess",
                "Parameters" : {
                    "model_part_name"         : "test.submodelpart_1",
                    "input_variable_name"     : "VELOCITY",
                    "output_variable_name"    : "VELOCITY",
                    "alignment_variable_name" : "NORMAL",
                    "is_tangential_alignment" : true,
                    "echo_level"              : 0
                }
            }
        ]''')

        node_ids = [1, 2, 3]
        check_values = []
        for node in self.model_part.Nodes:
            if (node.Id in node_ids):
                normal = node.GetSolutionStepValue(Kratos.NORMAL)
                normal_magnitude = float(
                    math.pow(
                        normal[0] * normal[0] + normal[1] * normal[1] +
                        normal[2] * normal[2], 0.5))
                unit_normal = [
                    normal[0] / normal_magnitude, normal[1] / normal_magnitude,
                    normal[2] / normal_magnitude
                ]
                vector = node.GetSolutionStepValue(Kratos.VELOCITY)
                vector_normal_magnitude = float(vector[0] * unit_normal[0] +
                                                vector[1] * unit_normal[1] +
                                                vector[2] * unit_normal[2])
                check_values.append([
                    vector[0] - vector_normal_magnitude * unit_normal[0],
                    vector[1] - vector_normal_magnitude * unit_normal[1],
                    vector[2] - vector_normal_magnitude * unit_normal[2]
                ])
            else:
                vector = node.GetSolutionStepValue(Kratos.VELOCITY)
                check_values.append(
                    [float(vector[0]),
                     float(vector[1]),
                     float(vector[2])])

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node, value in zip(self.model_part.Nodes, check_values):
            node_value = node.GetSolutionStepValue(Kratos.VELOCITY)
            for i in range(3):
                self.assertAlmostEqual(node_value[i], value[i], 9)

    def testVectorAlignProcessNormal(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyVectorAlignProcess",
                "Parameters" : {
                    "model_part_name"         : "test.submodelpart_1",
                    "input_variable_name"     : "VELOCITY",
                    "output_variable_name"    : "VELOCITY",
                    "alignment_variable_name" : "NORMAL",
                    "is_tangential_alignment" : false,
                    "echo_level"              : 0
                }
            }
        ]''')

        node_ids = [1, 2, 3]
        check_values = []
        for node in self.model_part.Nodes:
            if (node.Id in node_ids):
                normal = node.GetSolutionStepValue(Kratos.NORMAL)
                normal_magnitude = float(
                    math.pow(
                        normal[0] * normal[0] + normal[1] * normal[1] +
                        normal[2] * normal[2], 0.5))
                unit_normal = [
                    normal[0] / normal_magnitude, normal[1] / normal_magnitude,
                    normal[2] / normal_magnitude
                ]
                vector = node.GetSolutionStepValue(Kratos.VELOCITY)
                vector_normal_magnitude = float(vector[0] * unit_normal[0] +
                                                vector[1] * unit_normal[1] +
                                                vector[2] * unit_normal[2])
                check_values.append([
                    vector_normal_magnitude * unit_normal[0],
                    vector_normal_magnitude * unit_normal[1],
                    vector_normal_magnitude * unit_normal[2]
                ])
            else:
                vector = node.GetSolutionStepValue(Kratos.VELOCITY)
                check_values.append(
                    [float(vector[0]),
                     float(vector[1]),
                     float(vector[2])])

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node, value in zip(self.model_part.Nodes, check_values):
            node_value = node.GetSolutionStepValue(Kratos.VELOCITY)
            for i in range(3):
                self.assertAlmostEqual(node_value[i], value[i], 9)

    def testWallDistanceCalculationProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.submodelpart_1",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "STRUCTURE",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : "all"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSModellingApplication",
                "python_module" : "apply_custom_process",
                "process_name"  : "WallDistanceCalculationProcess",
                "Parameters" :             {
                    "model_part_name"          : "test",
                    "max_iterations"           : 1000,
                    "echo_level"               : 0,
                    "wall_flag_variable_name"  : "STRUCTURE",
                    "wall_flag_variable_value" : true,
                    "linear_solver_settings" : {
                        "solver_type"     : "amgcl"
                    }
                }
            }
        ]''')

        node_ids = [1, 2, 3, 4, 5, 6]
        check_values = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node, value in zip(self.model_part.Nodes, check_values):
            node_value = node.GetSolutionStepValue(Kratos.DISTANCE)
            self.assertAlmostEqual(node_value, value, 9)

# test model part is as follows
#       6   5   4
#       .---.---.
#      / |4 /\ 3|
#     / 1| /  \ |
#    /   |/  2 \|
#   .----.------.
#   1    2      3

    def __CreateModel(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        self.model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        self.model_part.CreateNewNode(5, 1.5, 1.0, 0.0)
        self.model_part.CreateNewNode(6, 1.0, 1.0, 0.0)
        prop = self.model_part.GetProperties()[0]

        self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 6], prop)
        self.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 5], prop)
        self.model_part.CreateNewElement("Element2D3N", 3, [3, 4, 5], prop)
        self.model_part.CreateNewElement("Element2D3N", 4, [2, 5, 6], prop)

        self.submodelpart_1 = self.model_part.CreateSubModelPart(
            "submodelpart_1")
        self.submodelpart_1.AddNodes([1, 2, 3])
        self.submodelpart_1.CreateNewCondition("Condition2D2N", 1, [1, 2],
                                               prop)
        self.submodelpart_1.CreateNewCondition("Condition2D2N", 2, [2, 3],
                                               prop)

        self.submodelpart_2 = self.model_part.CreateSubModelPart(
            "submodelpart_2")
        self.submodelpart_2.AddNodes([3, 4, 5])
        self.submodelpart_2.CreateNewCondition("Condition2D2N", 3, [3, 4],
                                               prop)
        self.submodelpart_2.CreateNewCondition("Condition2D2N", 4, [4, 5],
                                               prop)

        self.submodelpart_3 = self.model_part.CreateSubModelPart(
            "submodelpart_3")
        self.submodelpart_3.AddNodes([1, 6, 5])
        self.submodelpart_3.CreateNewCondition("Condition2D2N", 5, [5, 6],
                                               prop)
        self.submodelpart_3.CreateNewCondition("Condition2D2N", 6, [6, 1],
                                               prop)

        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        for node in self.model_part.Nodes:
            vector = Kratos.Vector(3)
            vector[0] = random.random()
            vector[1] = random.random()
            vector[2] = random.random()
            vector = node.SetSolutionStepValue(Kratos.VELOCITY, 0, vector)

            vector = Kratos.Vector(3)
            vector[0] = random.random()
            vector[1] = random.random()
            vector[2] = random.random()
            vector = node.SetSolutionStepValue(Kratos.NORMAL, 0, vector)

            scalar = random.random()
            node.SetSolutionStepValue(Kratos.DENSITY, 0, scalar)

    def __ExecuteProcesses(self):
        for process in self.process_list:
            process.Check()
        for process in self.process_list:
            process.ExecuteInitialize()
        for process in self.process_list:
            process.Execute()


if __name__ == '__main__':
    UnitTest.main()
