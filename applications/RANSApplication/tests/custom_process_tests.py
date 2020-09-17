import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities  import ReadModelPart


class CustomProcessTest(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        # add required variables to solution step list
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)

        cls.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
        cls.model_part.ProcessInfo.SetValue(Kratos.STEP, 1)

        ReadModelPart("BackwardFacingStepTest/backward_facing_step", cls.model_part)

    def setUp(self):
        # reinitialize variables for each test
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, Kratos.DENSITY, 0.0, 100.0, 0)
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, Kratos.VELOCITY, 0.0, 100.0, 0)
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, Kratos.PRESSURE, 0.0, 100.0, 0)
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, KratosRANS.TURBULENT_KINETIC_ENERGY, 0.0, 100.0, 0)
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, 0.0, 100.0, 0)
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, Kratos.KINEMATIC_VISCOSITY, 0.0, 100.0, 0)
        KratosRANS.RansTestUtilities.RandomFillNodalHistoricalVariable(self.model_part, Kratos.DISTANCE, 0.0, 100.0, 0)

    def testCheckScalarBoundsProcess(self):
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

    def testClipScalarVariableProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ClipScalarVariableProcess",
                "Parameters" : {
                    "model_part_name"                : "test",
                    "variable_name"                  : "DENSITY",
                    "min_value"                      : 20.0,
                    "max_value"                      : 60.0
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node in self.model_part.Nodes:
            density = node.GetSolutionStepValue(Kratos.DENSITY)
            self.assertEqual(density >= 20.0, True)
            self.assertEqual(density <= 60.0, True)

    def testApplyFlagProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.Slip2D.Slip2D_walls",
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
                    "model_part_name"                : "test.AutomaticInlet2D_inlet",
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

        for node in self.model.GetModelPart("test.AutomaticInlet2D_inlet").Nodes:
            self.assertEqual(node.Is(Kratos.INLET), True)

        for condition in self.model.GetModelPart("test.AutomaticInlet2D_inlet").Conditions:
            self.assertEqual(condition.Is(Kratos.INLET), True)
            self.assertEqual(condition.Is(Kratos.STRUCTURE), False)

        for node in self.model.GetModelPart("test.Slip2D.Slip2D_walls").Nodes:
            self.assertEqual(node.Is(Kratos.STRUCTURE), True)

        for condition in self.model.GetModelPart("test.Slip2D.Slip2D_walls").Conditions:
            self.assertEqual(condition.Is(Kratos.STRUCTURE), True)
            self.assertEqual(condition.Is(Kratos.INLET), False)

    def testLineOutputProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "LineOutputProcess",
                "Parameters" : {
                    "model_part_name"                   : "test",
                    "variable_names_list"               : ["DENSITY", "VELOCITY"],
                    "historical_value"                  : true,
                    "start_point"                       : [-0.09, 0.01, 0.0],
                    "end_point"                         : [0.19, 0.01, 0.0],
                    "number_of_sampling_points"         : 100,
                    "output_file_name"                  : "process_tests_data/line_output_test_output_historical",
                    "output_step_control_variable_name" : "STEP",
                    "output_step_interval"              : 1,
                    "write_header_information"          : false
                }
            },
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "LineOutputProcess",
                "Parameters" : {
                    "model_part_name"                   : "test",
                    "variable_names_list"               : [
                        "DENSITY",
                        "VELOCITY",
                        "EXTERNAL_FORCES_VECTOR",
                        "GREEN_LAGRANGE_STRAIN_TENSOR"],
                    "historical_value"                  : false,
                    "start_point"                       : [-0.09, 0.01, 0.0],
                    "end_point"                         : [0.19, 0.01, 0.0],
                    "number_of_sampling_points"         : 100,
                    "output_file_name"                  : "process_tests_data/line_output_test_output_non_historical",
                    "output_step_control_variable_name" : "STEP",
                    "output_step_interval"              : 1,
                    "write_header_information"          : false
                }
            },
            {
                "kratos_module" : "KratosMultiphysics",
                "python_module" : "compare_two_files_check_process",
                "Parameters" : {
                    "reference_file_name"   : "process_tests_data/line_output_test_output_historical_ref.csv",
                    "output_file_name"      : "process_tests_data/line_output_test_output_historical_1.000000.csv",
                    "remove_output_file"    : true,
                    "comparison_type"       : "csv_file",
                    "tolerance"             : 1e-6,
                    "relative_tolerance"    : 1e-9,
                    "dimension"             : 3
                }
            },
            {
                "kratos_module" : "KratosMultiphysics",
                "python_module" : "compare_two_files_check_process",
                "Parameters" : {
                    "reference_file_name"   : "process_tests_data/line_output_test_output_non_historical_ref.csv",
                    "output_file_name"      : "process_tests_data/line_output_test_output_non_historical_1.000000.csv",
                    "remove_output_file"    : true,
                    "comparison_type"       : "csv_file",
                    "tolerance"             : 1e-6,
                    "relative_tolerance"    : 1e-9,
                    "dimension"             : 3
                }
            }
        ]''')

        KratosRANS.RansTestUtilities.RandomFillNodalNonHistoricalVariable(self.model_part, Kratos.DENSITY, 0.0, 50.0)
        KratosRANS.RansTestUtilities.RandomFillNodalNonHistoricalVariable(self.model_part, Kratos.VELOCITY, 0.0, 50.0)

        for node in self.model_part.Nodes:
            v = Kratos.Vector(4)
            v[0] = node.X
            v[1] = node.Y
            v[2] = node.GetValue(Kratos.DENSITY)
            v[3] = node.GetValue(Kratos.DENSITY) * 1.1
            node.SetValue(Kratos.EXTERNAL_FORCES_VECTOR, v)

            m = Kratos.Matrix(2, 2)
            m[0, 0] = node.GetValue(Kratos.DENSITY)
            m[0, 1] = node.GetValue(Kratos.VELOCITY)[0]
            m[1, 0] = node.GetValue(Kratos.VELOCITY)[1]
            m[1, 1] = node.GetValue(Kratos.VELOCITY)[2]
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testKTurbulentIntensityInletProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "KTurbulentIntensityInletProcess",
                "Parameters" : {
                    "model_part_name"     : "test.AutomaticInlet2D_inlet",
                    "turbulent_intensity" : 0.01
                }
            }
        ]''')

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.model_part)

        test_variables = ["TURBULENT_KINETIC_ENERGY"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "k_turbulent_intensity_inlet_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node in self.model.GetModelPart(test_model_part_name).Nodes:
            self.assertEqual(node.IsFixed(KratosRANS.TURBULENT_KINETIC_ENERGY), True)

    def testEpsilonTurbulentMixingLengthInletProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "EpsilonTurbulentMixingLengthInletProcess",
                "Parameters" : {
                    "model_part_name"         : "test.AutomaticInlet2D_inlet",
                    "turbulent_mixing_length" : 0.005
                }
            }
        ]''')

        test_variables = ["TURBULENT_ENERGY_DISSIPATION_RATE"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "epsilon_turbulent_mixing_length_inlet_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.model_part)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node in self.model.GetModelPart(test_model_part_name).Nodes:
            self.assertEqual(node.IsFixed(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE), True)

    def testOmegaTurbulentMixingLengthInletProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "OmegaTurbulentMixingLengthInletProcess",
                "Parameters" : {
                    "model_part_name"     : "test"
                }
            }
        ]''')

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, self.model_part)

        test_variables = ["TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "omega_turbulent_mixing_length_inlet_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

        for node in self.model.GetModelPart(test_model_part_name).Nodes:
            self.assertEqual(node.IsFixed(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE), True)

    def testNutKEpsilonUpdateProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutKEpsilonUpdateProcess",
                "Parameters" : {
                    "model_part_name"     : "test"
                }
            }
        ]''')

        test_variables = ["TURBULENT_VISCOSITY", "VISCOSITY"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "nut_k_epsilon_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testNutKOmegaUpdateProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutKOmegaUpdateProcess",
                "Parameters" : {
                    "model_part_name"     : "test"
                }
            }
        ]''')

        test_variables = ["TURBULENT_VISCOSITY", "VISCOSITY"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "nut_k_omega_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testNutKOmegaSSTUpdateProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutKOmegaSSTUpdateProcess",
                "Parameters" : {
                    "model_part_name"     : "test"
                }
            }
        ]''')

        test_variables = ["TURBULENT_VISCOSITY", "VISCOSITY"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "nut_k_omega_sst_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testNutYPlusWallFunctionUpdateProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "NutYPlusWallFunctionUpdateProcess",
                "Parameters" : {
                    "model_part_name"     : "test"
                }
            }
        ]''')

        KratosRANS.RansTestUtilities.RandomFillConditionVariable(self.model.GetModelPart("test"), KratosRANS.RANS_Y_PLUS, 10.0, 100.0)

        test_variables = ["TURBULENT_VISCOSITY", "VISCOSITY"]
        test_model_part_name = "test.AutomaticInlet2D_inlet"
        test_file_name = "nut_y_plus_wall_function_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testComputeReactionsProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ComputeReactionsProcess",
                "Parameters" : {
                    "model_part_name"     : "test"
                }
            }
        ]''')

        KratosRANS.RansTestUtilities.RandomFillConditionVariable(self.model.GetModelPart("test"), KratosRANS.FRICTION_VELOCITY, 10.0, 100.0)

        test_variables = ["REACTION"]
        test_model_part_name = "test.Slip2D.Slip2D_walls"
        test_file_name = "compute_reactions_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def testWallDistanceCalculationProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ApplyFlagProcess",
                "Parameters" : {
                    "model_part_name"                : "test.Slip2D.Slip2D_walls",
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
                "Parameters" : {
                    "model_part_name"                  : "test",
                    "max_iterations"                   : 10,
                    "echo_level"                       : 0,
                    "wall_flag_variable_name"          : "STRUCTURE",
                    "wall_flag_variable_value"         : true,
                    "re_calculate_at_each_time_step"   : false,
                    "correct_distances_using_neighbors": true,
                    "linear_solver_settings" : {
                        "solver_type"     : "amgcl"
                    }
                }
            }
        ]''')

        test_variables = ["DISTANCE"]
        test_model_part_name = "test"
        test_file_name = "wall_distance_calculation_test_output"
        CustomProcessTest.__AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest.__AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__ExecuteProcesses()

    def __ExecuteProcesses(self):
        for process in self.process_list:
            process.Check()
        for process in self.process_list:
            process.ExecuteInitialize()
        for process in self.process_list:
            process.ExecuteBeforeSolutionLoop()
        for process in self.process_list:
            process.ExecuteInitializeSolutionStep()
        for process in self.process_list:
            if (hasattr(process, "ExecuteBeforeCouplingSolveStep")):
                process.ExecuteBeforeCouplingSolveStep()
        for process in self.process_list:
            if (hasattr(process, "ExecuteAfterCouplingSolveStep")):
                process.ExecuteAfterCouplingSolveStep()
        for process in self.process_list:
            process.ExecuteFinalizeSolutionStep()
        for process in self.process_list:
            process.ExecuteFinalize()

    @staticmethod
    def __AddJsonOutputProcess(settings, output_variables, output_model_part_name, output_file_name):
        settings_str = r"""
            {
                "kratos_module": "KratosMultiphysics",
                "python_module": "json_output_process",
                "process_name": "JsonOutputProcess",
                "Parameters": {
                    "output_variables": [
                        <VARIABLES_LIST>
                    ],
                    "output_file_name": "process_tests_data/<OUTPUT_FILE_NAME>.json",
                    "model_part_name": "<OUTPUT_MODEL_PART_NAME>",
                    "time_frequency": -2
                }
            }
        """
        settings_str = settings_str.replace("<VARIABLES_LIST>", CustomProcessTest.__GetVariablesString(output_variables))
        settings_str = settings_str.replace("<OUTPUT_FILE_NAME>", output_file_name)
        settings_str = settings_str.replace("<OUTPUT_MODEL_PART_NAME>", output_model_part_name)
        settings.Append(Kratos.Parameters(settings_str))

    @staticmethod
    def __AddJsonCheckProcess(settings, check_variables, model_part_name, input_file_name):
        settings_str = r"""
            {
                "kratos_module": "KratosMultiphysics",
                "python_module": "from_json_check_result_process",
                "help": "",
                "process_name": "FromJsonCheckResultProcess",
                "Parameters": {
                    "check_variables": [
                        <VARIABLES_LIST>
                    ],
                    "input_file_name": "process_tests_data/<INPUT_FILE_NAME>.json",
                    "model_part_name": "<MODEL_PART_NAME>",
                    "tolerance": 1e-9,
                    "relative_tolerance": 1e-12,
                    "time_frequency": -2
                }
            }
            """
        settings_str = settings_str.replace("<VARIABLES_LIST>", CustomProcessTest.__GetVariablesString(check_variables))
        settings_str = settings_str.replace("<INPUT_FILE_NAME>", input_file_name + "_ref")
        settings_str = settings_str.replace("<MODEL_PART_NAME>", model_part_name)
        settings.Append(Kratos.Parameters(settings_str))

    @staticmethod
    def __GetVariablesString(variables_list):
        for i, variable in enumerate(variables_list):
            variables_list[i] = "\"{:s}\"".format(variable)
        return ",".join(variables_list)

if __name__ == '__main__':
    UnitTest.main()
