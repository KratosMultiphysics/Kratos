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

    def testLineOutputProcess(self):
        self.__CreateModel()
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "LineOutputProcess",
                "Parameters" : {
                    "model_part_name"                : "test",
                    "variable_names_list"               : ["DENSITY", "VELOCITY"],
                    "historical_value"                  : true,
                    "start_point"                       : [0.0, 0.0, 0.0],
                    "end_point"                         : [1.0, 0.0, 0.0],
                    "number_of_sampling_points"         : 100,
                    "output_file_name"                  : "process_tests_data/line_output_test_output",
                    "output_step_control_variable_name" : "STEP",
                    "output_step_interval"              : 1,
                    "write_header_information"          : false
                }
            },
            {
                "kratos_module" : "KratosMultiphysics",
                "python_module" : "compare_two_files_check_process",
                "Parameters" : {
                    "reference_file_name"   : "process_tests_data/line_output_test_output_ref.csv",
                    "output_file_name"      : "process_tests_data/line_output_test_output_2.csv",
                    "remove_output_file"    : true,
                    "comparison_type"       : "csv_file",
                    "tolerance"             : 1e-6,
                    "relative_tolerance"    : 1e-9,
                    "dimension"             : 3
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()


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
        self.model_part.CloneTimeStep(1.0)
        self.model_part.ProcessInfo[Kratos.STEP] = 2

        random.seed("rans_test")
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
