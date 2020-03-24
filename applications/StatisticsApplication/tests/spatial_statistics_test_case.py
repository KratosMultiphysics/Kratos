from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
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

    def testSpatialStatisticsProcessNodalHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultSettings("nodal_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessNodalNonHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultSettings("nodal_non_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessElementNonHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultSettings("element_non_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessConditionNonHistorical(self):
        settings = SpatialStatisticsTestCase.__GetDefaultSettings("condition_non_historical")
        self.__TestMethod(settings)

    def __TestMethod(self, settings):
        factory = KratosProcessFactory(self.current_model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        InitializeProcesses(self)

        for step in range(0, 12, 2):
            self.model_part.CloneTimeStep(step)
            self.model_part.ProcessInfo[Kratos.STEP] = step
            InitializeModelPartVariables(self.model_part, False)
            ExecuteProcessFinalizeSolutionStep(self)

        for process in self.process_list:
            process.ExecuteFinalize()

    @staticmethod
    def __GetDefaultSettings(container_name):
        settings_str = r'''
            {
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "spatial_statistics_process",
                "Parameters" : {
                    "model_part_name" : "test_model_part",
                    "input_variable_settings" : [
                        {
                            "method_name"    : "sum",
                            "norm_type"      : "none",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "mean",
                            "norm_type"      : "none",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "variance",
                            "norm_type"      : "none",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "rootmeansquare",
                            "norm_type"      : "none",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "sum",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "mean",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "variance",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "rootmeansquare",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "min",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "max",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "median",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "distribution",
                            "norm_type"      : "magnitude",
                            "container"      : "<CONTAINER_NAME>",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        }
                    ],
                    "output_settings" : {
                        "output_control_variable": "STEP",
                        "output_time_interval"   : 1,
                        "write_kratos_version"   : false,
                        "write_time_stamp"       : false,
                        "output_file_settings"   : {
                            "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
                            "folder_name": "spatial_statistics_process",
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
                    "output_file_name": "spatial_statistics_process/<OUTPUT_FILE_NAME>",
                    "reference_file_name": "spatial_statistics_process/<INPUT_FILE_NAME>",
                    "comparison_type": "dat_file",
                    "remove_output_file": true,
                    "tolerance": 1e-16
                }
            }
        '''
        final_settings_str = "[\n" + settings_str.replace("<CONTAINER_NAME>", container_name)
        method_list = {}
        method_list["sum"] = ["none", "magnitude"]
        method_list["mean"] = ["none", "magnitude"]
        method_list["variance"] = ["none", "magnitude"]
        method_list["rootmeansquare"] = ["none", "magnitude"]
        method_list["min"] = ["magnitude"]
        method_list["max"] = ["magnitude"]
        method_list["median"] = ["magnitude"]
        method_list["distribution"] = ["magnitude"]

        for method, norm_types in method_list.items():
            for norm_type in norm_types:
                name_prefix = "test_model_part_" + container_name + "_" + norm_type + "_" + method
                file_output_name = name_prefix + ".dat"
                reference_file_name = name_prefix + ".ref.dat"
                current_settings = check_settings.replace("<OUTPUT_FILE_NAME>", file_output_name)
                current_settings = current_settings.replace("<INPUT_FILE_NAME>", reference_file_name)
                final_settings_str += ",\n" + current_settings

        final_settings_str += "\n]"
        return Kratos.Parameters(final_settings_str)
