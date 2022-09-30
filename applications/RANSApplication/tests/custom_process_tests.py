from math import sqrt

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities  import ReadModelPart

from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.FluidDynamicsApplication.check_and_prepare_model_process_fluid import CheckAndPrepareModelProcess


class CustomProcessTest(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("FluidModelPart")
        cls.model_part.SetBufferSize(2)

        # add required variables to solution step list
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
        cls.model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)

        cls.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
        cls.model_part.ProcessInfo.SetValue(Kratos.STEP, 1)
        cls.model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.1)

        with UnitTest.WorkFolderScope(".", __file__):
            ReadModelPart("BackwardFacingStepTest/backward_facing_step", cls.model_part)
            CheckAndPrepareModelProcess(cls.model_part,
                                        Kratos.Parameters("""{
                "volume_model_part_name": "Parts_fluid",
                "skin_parts" : ["AutomaticInlet2D_inlet", "Outlet2D_outlet", "Slip2D"],
                "assign_neighbour_elements_to_conditions": true
            }""")).Execute()

    def setUp(self):
        # reinitialize variables for each test
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, Kratos.VELOCITY, 0.0, 100.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, Kratos.ACCELERATION, 0.0, 10.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, Kratos.ACCELERATION, 20.0, 50.0, 1)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, Kratos.PRESSURE, 0.0, 100.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, KratosRANS.TURBULENT_KINETIC_ENERGY, 0.0, 100.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 20.0, 100.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 50.0, 100.0, 1)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, 10.0, 100.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillHistoricalVariable(self.model_part, Kratos.DISTANCE, 0.0, 100.0, 0)

        Kratos.VariableUtils().SetVariable(Kratos.DENSITY, 1.0, self.model_part.Nodes)

    def testCheckScalarBoundsProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "CheckScalarBoundsProcess",
                "Parameters" : {
                    "model_part_name"                : "FluidModelPart",
                    "variable_name"                  : "DENSITY"
                }
            }
        ]''')

        self._RunProcessTest(settings)

    def testClipScalarVariableProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "ClipScalarVariableProcess",
                "Parameters" : {
                    "model_part_name"                : "FluidModelPart",
                    "variable_name"                  : "DENSITY",
                    "min_value"                      : 20.0,
                    "max_value"                      : 60.0
                }
            }
        ]''')

        self._RunProcessTest(settings)

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
                    "model_part_name"                : "FluidModelPart.Slip2D.Slip2D_walls",
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
                    "model_part_name"                : "FluidModelPart.AutomaticInlet2D_inlet",
                    "echo_level"                     : 0,
                    "flag_variable_name"             : "INLET",
                    "flag_variable_value"            : true,
                    "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
                }
            }
        ]''')

        self._RunProcessTest(settings)

        for node in self.model.GetModelPart("FluidModelPart.AutomaticInlet2D_inlet").Nodes:
            self.assertEqual(node.Is(Kratos.INLET), True)

        for condition in self.model.GetModelPart("FluidModelPart.AutomaticInlet2D_inlet").Conditions:
            self.assertEqual(condition.Is(Kratos.INLET), True)
            self.assertEqual(condition.Is(Kratos.STRUCTURE), False)

        for node in self.model.GetModelPart("FluidModelPart.Slip2D.Slip2D_walls").Nodes:
            self.assertEqual(node.Is(Kratos.STRUCTURE), True)

        for condition in self.model.GetModelPart("FluidModelPart.Slip2D.Slip2D_walls").Conditions:
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
                    "model_part_name"                   : "FluidModelPart",
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
                    "model_part_name"                   : "FluidModelPart",
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

        KratosCFD.FluidTestUtilities.RandomFillNonHistoricalVariable(self.model_part.Nodes, Kratos.DENSITY, 2, 0.0, 50.0)
        KratosCFD.FluidTestUtilities.RandomFillNonHistoricalVariable(self.model_part.Nodes, Kratos.VELOCITY, 2, 0.0, 50.0)

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

        self._RunProcessTest(settings)

    def testKTurbulentIntensityInletProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "KTurbulentIntensityInletProcess",
                "Parameters" : {
                    "model_part_name"     : "FluidModelPart.AutomaticInlet2D_inlet",
                    "turbulent_intensity" : 0.01
                }
            }
        ]''')

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.model_part)

        test_variables = ["TURBULENT_KINETIC_ENERGY"]
        test_model_part_name = "FluidModelPart.AutomaticInlet2D_inlet"
        test_file_name = "k_turbulent_intensity_inlet_test_output"
        CustomProcessTest._AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest._AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        self._RunProcessTest(settings)

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
                    "model_part_name"         : "FluidModelPart.AutomaticInlet2D_inlet",
                    "turbulent_mixing_length" : 0.005
                }
            }
        ]''')

        test_variables = ["TURBULENT_ENERGY_DISSIPATION_RATE"]
        test_model_part_name = "FluidModelPart.AutomaticInlet2D_inlet"
        test_file_name = "epsilon_turbulent_mixing_length_inlet_test_output"
        CustomProcessTest._AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest._AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.model_part)
        self.model_part.ProcessInfo.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, 0.09)

        self._RunProcessTest(settings)

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
                    "model_part_name"     : "FluidModelPart"
                }
            }
        ]''')

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, self.model_part)
        self.model_part.ProcessInfo.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, 0.09)

        test_variables = ["TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE"]
        test_model_part_name = "FluidModelPart.AutomaticInlet2D_inlet"
        test_file_name = "omega_turbulent_mixing_length_inlet_test_output"
        CustomProcessTest._AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest._AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        self._RunProcessTest(settings)

        for node in self.model.GetModelPart(test_model_part_name).Nodes:
            self.assertEqual(node.IsFixed(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE), True)

    def testWallDistanceCalculationProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "WallDistanceCalculationProcess",
                "Parameters" : {
                    "main_model_part_name"                 : "FluidModelPart",
                    "wall_model_part_name"                 : "FluidModelPart.Slip2D.Slip2D_walls",
                    "echo_level"                           : 0,
                    "max_distance"                         : 1e+30,
                    "max_levels"                           : 14,
                    "wall_distance_update_execution_points": ["initialize"]
                }
            }
        ]''')

        test_variables = ["DISTANCE"]
        test_model_part_name = "FluidModelPart"
        test_file_name = "wall_distance_calculation_test_output"
        CustomProcessTest._AddJsonCheckProcess(settings, test_variables, test_model_part_name, test_file_name)
        # CustomProcessTest._AddJsonOutputProcess(settings, test_variables, test_model_part_name, test_file_name)

        CalculateNormalsOnConditions(self.model_part)
        self._RunProcessTest(settings)

    def testVariableDataTransferProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "VariableDataTransferProcess",
                "Parameters" : {
                    "source_model_part_name"     : "FluidModelPart",
                    "destination_model_part_name": "FluidModelPart",
                    "copy_execution_points"      : ["initialize"],
                    "copy_variables_list"        : [
                        {
                            "source_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "ACCELERATION",
                                "step_index"             : 0
                            },
                            "destination_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "REACTION",
                                "step_index"             : 0
                            }
                        },
                        {
                            "source_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "ACCELERATION",
                                "step_index"             : 1
                            },
                            "destination_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "REACTION",
                                "step_index"             : 1
                            }
                        },
                        {
                            "source_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "ACCELERATION",
                                "step_index"             : 0
                            },
                            "destination_variable_settings" : {
                                "is_historical_container": false,
                                "variable_name"          : "ACCELERATION",
                                "step_index"             : 0
                            }
                        },
                        {
                            "source_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "ACCELERATION",
                                "step_index"             : 1
                            },
                            "destination_variable_settings" : {
                                "is_historical_container": false,
                                "variable_name"          : "NORMAL",
                                "step_index"             : 0
                            }
                        },
                        {
                            "source_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "REACTION",
                                "step_index"             : 0
                            },
                            "destination_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "VELOCITY",
                                "step_index"             : 1
                            }
                        },
                        {
                            "source_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "ACCELERATION",
                                "step_index"             : 1
                            },
                            "destination_variable_settings" : {
                                "is_historical_container": true,
                                "variable_name"          : "VELOCITY",
                                "step_index"             : 0
                            }
                        }
                    ],
                    "echo_level"                 : 0
                }
            }
        ]''')

        self._RunProcessTest(settings)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 1), node.GetSolutionStepValue(Kratos.REACTION, 1), 9)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 0), node.GetSolutionStepValue(Kratos.REACTION, 0), 9)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 1), node.GetValue(Kratos.NORMAL), 9)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 0), node.GetValue(Kratos.ACCELERATION), 9)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 1), node.GetSolutionStepValue(Kratos.VELOCITY, 0), 9)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 0), node.GetSolutionStepValue(Kratos.VELOCITY, 1), 9)

    def testInitializeBossakPreviousStepVariableDerivatives(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "InitializeBossakPreviousStepVariableDerivatives",
                "Parameters" : {
                    "model_part_name" : "FluidModelPart",
                    "echo_level"      : 1,
                    "variables_list"  : [
                        {
                            "first_derivative_variable_name"         : "TURBULENT_KINETIC_ENERGY",
                            "second_derivative_variable_name"        : "TURBULENT_KINETIC_ENERGY_RATE",
                            "relaxed_second_derivative_variable_name": "RANS_AUXILIARY_VARIABLE_1",
                            "is_relaxed_second_derivative_historical": true
                        }
                    ]
                }
            }
        ]''')

        scalar_scheme = KratosRANS.BossakRelaxationScalarScheme(-0.3, 1.0, KratosRANS.TURBULENT_KINETIC_ENERGY)
        scalar_scheme.InitializeSolutionStep(self.model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())
        scalar_scheme.UpdateScalarRateVariables(self.model_part)

        second_derivative_values = {}
        for node in self.model_part.Nodes:
            second_derivative_values[node.Id] = [
                node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 0),
                node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 1)]
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 1, 0.0)

        self._RunProcessTest(settings)

        for node in self.model_part.Nodes:
            current_nodal_data = second_derivative_values[node.Id]
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 0), current_nodal_data[0], 9)
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, 1), current_nodal_data[1], 9)

    def testInitializePreviousSolutionStepValuesProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "initialize_previous_solution_step_values_process",
                "process_name"  : "InitializePreviousSolutionStepValuesProcess",
                "Parameters" : {
                    "model_part_name" : "FluidModelPart",
                    "echo_level"      : 0,
                    "variable_name"   : "ACCELERATION_X",
                    "value"           : "x+y*t"
                }
            }
        ]''')

        self.model_part.CloneTimeStep(0.5)
        self.model_part.CloneTimeStep(2.0)

        current_values = {}
        for node in self.model_part.Nodes:
            current_values[node.Id] = [node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.ACCELERATION, 1)]

        self._GetProcessList(settings)
        self.process_list[0].Execute()

        for node in self.model_part.Nodes:
            data = current_values[node.Id]
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), data[0], 9)
            v =  data[1]
            v[0] = node.X + node.Y * 0.5
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION, 1), v, 9)

    def testCheckScalarBoundsProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.RANSApplication",
                "python_module" : "cpp_process_factory",
                "process_name"  : "SmoothClipScalarVariableProcess",
                "Parameters" : {
                    "model_part_name"                           : "FluidModelPart",
                    "variable_name"                             : "PRESSURE",
                    "echo_level"                                : 1,
                    "min_value"                                 : 1e-18,
                    "max_value"                                 : 300,
                    "max_number_of_sweeps"                      : 10000,
                    "inverse_distance_weighting_power_parameter": 2,
                    "always_find_neighbour_nodes"               : false,
                    "execution_points"                          : ["after_coupling_solve_step"]
                }
            }
        ]''')

        neighbour_finding_process = Kratos.FindGlobalNodalNeighboursProcess(self.model_part.GetCommunicator().GetDataCommunicator(), self.model_part)
        neighbour_finding_process.Execute()
        neighbour_id_map = neighbour_finding_process.GetNeighbourIds(self.model_part.Nodes)

        node_id_list = [229, 216, 221, 224, 233, 236, 241]
        for node_id in node_id_list:
            self.model_part.GetNode(node_id).SetSolutionStepValue(Kratos.PRESSURE, -10)

        reference_values = {}
        for node in self.model_part.Nodes:
            reference_values[node.Id] = node.GetSolutionStepValue(Kratos.PRESSURE)

        def smoothing_layer(node_ids_list):
            for smoothing_layer_node_id in node_ids_list:
                smoothing_layer_node = self.model_part.GetNode(smoothing_layer_node_id)
                total_weight = 0.0
                current_value = 0.0
                neighbour_nodes_list = neighbour_id_map[smoothing_layer_node_id]
                for neighbour_node_id in neighbour_nodes_list:
                    neighbour_node = self.model_part.GetNode(neighbour_node_id)
                    neighbour_node_value = neighbour_node.GetSolutionStepValue(Kratos.PRESSURE)
                    if neighbour_node_value > 1e-18 and neighbour_node_value < 300 and not neighbour_node_id in node_ids_list:
                        weight = (neighbour_node.X-smoothing_layer_node.X)**2 + (neighbour_node.Y-smoothing_layer_node.Y)**2 + (neighbour_node.Z-smoothing_layer_node.Z)**2
                        current_value += neighbour_node_value / weight
                        total_weight += 1/weight
                smoothing_layer_node.SetSolutionStepValue(Kratos.PRESSURE, current_value/total_weight)

        smoothing_layer([216, 221, 224, 233, 236, 241])
        smoothing_layer([229])

        for node_id in node_id_list:
            print(node_id)
            reference_values[node_id] = self.model_part.GetNode(node_id).GetSolutionStepValue(Kratos.PRESSURE)
            self.model_part.GetNode(node_id).SetSolutionStepValue(Kratos.PRESSURE, -10)

        self._RunProcessTest(settings)

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.PRESSURE), reference_values[node.Id], 9)

    def _GetProcessList(self, settings):
        with UnitTest.WorkFolderScope(".", __file__):
            factory = KratosProcessFactory(self.model)
            self.process_list = factory.ConstructListOfProcesses(settings)

    def _RunProcessTest(self, settings):
        self._GetProcessList(settings)
        with UnitTest.WorkFolderScope(".", __file__):
            self._ExecuteProcesses()

    def _ExecuteProcesses(self):
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
    def _AddJsonOutputProcess(settings, output_variables, output_model_part_name, output_file_name):
        settings_str = r"""
            {
                "kratos_module": "KratosMultiphysics",
                "python_module": "json_output_process",
                "process_name": "JsonOutputProcess",
                "Parameters": {
                    "output_variables": [
                        <VARIABLES_LIST>
                    ],
                    "output_file_name": "process_tests_data/<OUTPUT_FILE_NAME>_ref.json",
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
    def _AddJsonCheckProcess(settings, check_variables, model_part_name, input_file_name):
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
