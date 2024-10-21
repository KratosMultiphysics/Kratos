import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.process_factory import KratosProcessFactory
import KratosMultiphysics.StatisticsApplication as KratosStats

from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeModelPartVariables
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeProcesses
from KratosMultiphysics.StatisticsApplication.test_utilities import ExecuteProcessFinalizeSolutionStep

import statistics_test_case

class SpatialStatisticsTestCase(statistics_test_case.StatisticsTestCase):
    @classmethod
    def AddVariables(cls):
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        # output variables for output_1
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.SCALAR_NORM)
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.VECTOR_3D_NORM)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.YIELD_STRESS)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.CUTTED_AREA)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NET_INPUT_MATERIAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.WET_VOLUME)

    def setUp(self):
        self.model_part.ProcessInfo[Kratos.TIME] = 0.0
        self.model_part.ProcessInfo[Kratos.STEP] = 0

    def testSpatialStatisticsProcessNodalHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultParameters("nodal_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessNodalNonHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultParameters("nodal_non_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessElementNonHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultParameters("element_non_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessConditionNonHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultParameters("condition_non_historical")
        self.__TestMethod(settings)

    def __TestMethod(self, settings):
        with KratosUnittest.WorkFolderScope(".", __file__):
            factory = KratosProcessFactory(self.current_model)
            self.process_list = factory.ConstructListOfProcesses(settings)
            InitializeProcesses(self)

            for step in range(0, 12, 2):
                self.model_part.CloneTimeStep(step)
                self.model_part.ProcessInfo[Kratos.STEP] = step
                InitializeModelPartVariables(self.model_part, False)
                ExecuteProcessFinalizeSolutionStep(self)

                for process in self.process_list:
                    if isinstance(process, Kratos.OutputProcess) and process.IsOutputStep():
                        process.PrintOutput()

            for process in self.process_list:
                process.ExecuteFinalize()

    @staticmethod
    def __GetDefaultParameters(container_name):
        settings_str = r'''
            {
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "spatial_statistics_process",
                "Parameters" : {
                    "model_part_name" : "test_model_part",
                    "echo_level"      : 0,
                    "input_variable_settings" : [
                        {
                            "norm_type"      : "none",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY"]
                        },
                        {
                            "norm_type"      : "l2",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"]
                        }
                    ],
                    "statistics_methods": [
                        {
                            "method_name"    : "sum"
                        },
                        {
                            "method_name"    : "mean"
                        },
                        {
                            "method_name"    : "variance"
                        },
                        {
                            "method_name"    : "rootmeansquare"
                        },
                        {
                            "method_name"    : "min"
                        },
                        {
                            "method_name"    : "max"
                        },
                        {
                            "method_name"    : "median"
                        },
                        {
                            "method_name"    : "distribution",
                            "method_settings": {
                                "min_value": -1.0,
                                "max_value": 100.0
                            }
                        }
                    ],
                    "output_settings" : {
                        "output_control_variable": "STEP",
                        "output_time_interval"   : 2,
                        "write_kratos_version"   : false,
                        "write_time_stamp"       : false,
                        "output_file_settings"   : {
                            "file_name"  : "<model_part_name>_<CONTAINER_NAME>.dat",
                            "output_path": "spatial_statistics_process",
                            "write_buffer_size" : -1
                        }
                    }
                }
            }
        '''
        check_settings = r'''
            {
                "python_module": "compare_two_files_check_process",
                "kratos_module": "KratosMultiphysics",
                "help": "",
                "process_name": "CompareTwoFilesCheckProcess",
                "Parameters": {
                    "output_file_name": "spatial_statistics_process/test_model_part_<CONTAINER_NAME>.dat",
                    "reference_file_name": "spatial_statistics_process/test_model_part_<CONTAINER_NAME>.ref.dat",
                    "comparison_type": "deterministic",
                    "remove_output_file": true,
                    "tolerance": 1e-16
                }
            }
        '''
        final_settings_str = "[\n" + settings_str.replace("<CONTAINER_NAME>", container_name) + ",\n" + check_settings.replace("<CONTAINER_NAME>", container_name) + "\n]"
        return Kratos.Parameters(final_settings_str)
