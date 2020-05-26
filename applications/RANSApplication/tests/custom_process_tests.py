import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.KratosUnittest as UnitTest
import random, math


class CustomProcessTest(UnitTest.TestCase):
    def testCheckScalarBoundsProcess(self):
        self.__CreateModel()
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "CheckScalarBoundsProcess",
                "Parameters" : {
                    "model_part_name"                : "test",
                    "variable_name"                  : "DENSITY"
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testCheckVectorBoundsProcess(self):
        self.__CreateModel()
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "CheckVectorBoundsProcess",
                "Parameters" : {
                    "model_part_name"                : "test",
                    "variable_name"                  : "VELOCITY"
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testClipScalarVariableProcess(self):
        self.__CreateModel()
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ClipScalarVariableProcess",
                "Parameters" : {
                    "model_part_name"                : "test",
                    "variable_name"                  : "DENSITY",
                    "min_value"                      : 0.5,
                    "max_value"                      : 0.6
                }
            }
        ]''')

        model_part = self.model["test"]
        for node in model_part.Nodes:
            density = node.GetSolutionStepValue(Kratos.DENSITY)
            node.SetSolutionStepValue(Kratos.DENSITY, 0, math.pow(-1, node.Id) * density)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node in model_part.Nodes:
            density = node.GetSolutionStepValue(Kratos.DENSITY)
            self.assertEqual(density >= 0.5, True)
            self.assertEqual(density <= 0.6, True)

    def testApplyFlagProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.submodelpart_1",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "STRUCTURE",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.submodelpart_2",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "INLET",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
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
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "FindNodalNeighboursProcess",
                "Parameters" : {
                    "model_part_name"      : "test"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "FindConditionParentProcess",
                "Parameters" : {
                    "model_part_name"      : "test.submodelpart_1"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ScalarCellCenteredAveragingProcess",
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
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "FindNodalNeighboursProcess",
                "Parameters" : {
                    "model_part_name"      : "test"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "FindConditionParentProcess",
                "Parameters" : {
                    "model_part_name"      : "test.submodelpart_1"
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "VectorCellCenteredAveragingProcess",
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
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "VectorAlignProcess",
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
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "VectorAlignProcess",
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
        self.__CreateModel(
            variable_list=[Kratos.DISTANCE, Kratos.FLAG_VARIABLE])

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.submodelpart_1",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "STRUCTURE",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
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

        check_values = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node, value in zip(self.model_part.Nodes, check_values):
            node_value = node.GetSolutionStepValue(Kratos.DISTANCE)
            self.assertAlmostEqual(node_value, value, 9)

    def testLogarithmicYPlusCalculationProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "LogarithmicYPlusCalculationProcess",
                "Parameters" :             {
                    "model_part_name" : "test",
                    "echo_level"      : 0,
                    "max_iterations"  : 10,
                    "tolerance"       : 1e-6,
                    "constants": {
                        "von_karman"  : 0.41,
                        "beta"        : 5.2
                    }
                }
            }
        ]''')

        # for node, distance in zip(self.model_part.Nodes, distance_value):
        #     node.SetSolutionStepValue(Kratos.DISTANCE, 0, distance)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        beta = 5.2
        von_karman = 0.41
        for node in self.model_part.Nodes:
            y_plus = node.GetSolutionStepValue(KratosRANS.RANS_Y_PLUS)
            nu = node.GetSolutionStepValue(Kratos.KINEMATIC_VISCOSITY)
            y = node.GetSolutionStepValue(Kratos.DISTANCE)
            u = node.GetSolutionStepValue(Kratos.VELOCITY)
            if (y > 0.0):
                u_tau = y_plus * nu / y
                u_plus = math.sqrt(u[0] * u[0] + u[1] * u[1] +
                                   u[2] * u[2]) / u_tau
                if (y_plus > 11.06):
                    self.assertAlmostEqual(
                        u_plus, 1 / von_karman * math.log(y_plus) + beta, 9)
                    print("Checked logarithmic region...")
                else:
                    self.assertAlmostEqual(u_plus, y_plus, 9)
                    print("Checked linear region...")
            else:
                self.assertAlmostEqual(y_plus, 0.0, 9)

    def testLogarithmicYPlusVelocitySensitivitiesProcessFlow(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "LogarithmicYPlusCalculationProcess",
                "Parameters" :             {
                    "model_part_name" : "test",
                    "echo_level"      : 0,
                    "max_iterations"  : 10,
                    "tolerance"       : 1e-6,
                    "constants": {
                        "von_karman"  : 0.41,
                        "beta"        : 5.2
                    }
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "LogarithmicYPlusVelocitySensitivitiesProcess",
                "Parameters" :             {
                    "model_part_name" : "test",
                    "echo_level"      : 0,
                    "constants": {
                        "von_karman"  : 0.41,
                        "beta"        : 5.2
                    }
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testNutKEpsilonHighReCalculationProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutKEpsilonHighReCalculationProcess",
                "Parameters" :             {
                    "model_part_name" : "test",
                    "echo_level"      : 0,
                    "c_mu"            : 0.09
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        c_mu = 0.09
        for node in self.model_part.Nodes:
            k = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY)
            epsilon = node.GetSolutionStepValue(
                KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
            node_value = node.GetSolutionStepValue(Kratos.TURBULENT_VISCOSITY)
            self.assertAlmostEqual(node_value, c_mu * k * k / epsilon, 9)

    def testNutKEpsilonHighReSensitivitiesProcess(self):
        self.__CreateModel()

        adjoint_settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutKEpsilonHighReSensitivitiesProcess",
                "Parameters" :             {
                    "model_part_name" : "test"
                }
            }
        ]''')

        primal_settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutKEpsilonHighReCalculationProcess",
                "Parameters" :             {
                    "model_part_name" : "test"
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(adjoint_settings)
        self.__ExecuteProcesses()

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(primal_settings)
        self.__ExecuteProcesses()

        delta = 1e-9
        variable_list = [
            KratosRANS.TURBULENT_KINETIC_ENERGY,
            KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE
        ]
        check_variable = Kratos.TURBULENT_VISCOSITY
        sensitivity_variable = KratosRANS.RANS_NUT_SCALAR_PARTIAL_DERIVATIVES

        for node in self.model_part.Nodes:
            for (variable, index) in zip(variable_list,
                                         range(len(variable_list))):
                self.__ExecuteProcesses()
                current_check_value = node.GetSolutionStepValue(check_variable)

                current_value = node.GetSolutionStepValue(variable)
                node.SetSolutionStepValue(variable, 0, current_value + delta)

                self.__ExecuteProcesses()
                perturbed_check_value = node.GetSolutionStepValue(
                    check_variable)

                node.SetSolutionStepValue(variable, 0, current_value)

                check_value_sensitivity = (
                    perturbed_check_value - current_check_value) / delta
                adjoint_sensitivity = node.GetValue(
                    sensitivity_variable)[index]

                self.assertTrue(UnitTest.isclose(adjoint_sensitivity, check_value_sensitivity, rel_tol=1e-6))


# test model part is as follows
#       6   5   4
#       .---.---.
#      / |4 /\ 3|
#     / 1| /  \ |
#    /   |/  2 \|
#   .----.------.
#   1    2      3

    def __CreateModel(
            self,
            variable_list=[
                Kratos.VELOCITY, Kratos.NORMAL, Kratos.DENSITY,
                Kratos.KINEMATIC_VISCOSITY, Kratos.DISTANCE,
                Kratos.FLAG_VARIABLE, Kratos.TURBULENT_VISCOSITY,
                KratosRANS.RANS_Y_PLUS, KratosRANS.TURBULENT_KINETIC_ENERGY,
                KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE
            ]):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        for variable in variable_list:
            self.model_part.AddNodalSolutionStepVariable(variable)

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
        self.submodelpart_1.CreateNewCondition("LineCondition2D2N", 1, [1, 2],
                                               prop)
        self.submodelpart_1.CreateNewCondition("LineCondition2D2N", 2, [2, 3],
                                               prop)

        self.submodelpart_2 = self.model_part.CreateSubModelPart(
            "submodelpart_2")
        self.submodelpart_2.AddNodes([3, 4, 5])
        self.submodelpart_2.CreateNewCondition("LineCondition2D2N", 3, [3, 4],
                                               prop)
        self.submodelpart_2.CreateNewCondition("LineCondition2D2N", 4, [4, 5],
                                               prop)

        self.submodelpart_3 = self.model_part.CreateSubModelPart(
            "submodelpart_3")
        self.submodelpart_3.AddNodes([1, 6, 5])
        self.submodelpart_3.CreateNewCondition("LineCondition2D2N", 5, [5, 6],
                                               prop)
        self.submodelpart_3.CreateNewCondition("LineCondition2D2N", 6, [6, 1],
                                               prop)

        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        self.model_part.ProcessInfo[KratosRANS.TURBULENCE_RANS_C_MU] = 0.09

        for node in self.model_part.Nodes:
            vector = Kratos.Vector(3)
            vector[0] = random.random()
            vector[1] = random.random()
            vector[2] = random.random()
            if (node.SolutionStepsDataHas(Kratos.VELOCITY)):
                vector = node.SetSolutionStepValue(Kratos.VELOCITY, 0, vector)

            vector = Kratos.Vector(3)
            vector[0] = random.random()
            vector[1] = random.random()
            vector[2] = random.random()
            if (node.SolutionStepsDataHas(Kratos.NORMAL)):
                vector = node.SetSolutionStepValue(Kratos.NORMAL, 0, vector)

            scalar = random.random()
            if (node.SolutionStepsDataHas(Kratos.DENSITY)):
                node.SetSolutionStepValue(Kratos.DENSITY, 0, scalar)

            scalar = random.random() * 1e-3
            if (node.SolutionStepsDataHas(Kratos.KINEMATIC_VISCOSITY)):
                node.SetSolutionStepValue(Kratos.KINEMATIC_VISCOSITY, 0,
                                          scalar)

            scalar = random.random()
            if (node.SolutionStepsDataHas(
                    KratosRANS.TURBULENT_KINETIC_ENERGY)):
                node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY,
                                          0, scalar)
            scalar = random.random()
            if (node.SolutionStepsDataHas(
                    KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)):
                node.SetSolutionStepValue(
                    KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, 0, scalar)

            scalar = random.random()
            if (node.SolutionStepsDataHas(Kratos.DISTANCE)):
                node.SetSolutionStepValue(Kratos.DISTANCE, 0, scalar)

    def __ExecuteProcesses(self):
        for process in self.process_list:
            process.Check()
        for process in self.process_list:
            process.ExecuteInitialize()
        for process in self.process_list:
            process.Execute()

if __name__ == '__main__':
    UnitTest.main()
