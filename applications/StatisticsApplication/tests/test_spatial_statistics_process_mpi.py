import KratosMultiphysics as Kratos

import KratosMultiphysics as Kratos
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.mpi as KratosMPI
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.StatisticsApplication as KratosStats
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.StatisticsApplication.test_utilities import CreateModelPart
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeModelPartVariables
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeProcesses
from KratosMultiphysics.StatisticsApplication.test_utilities import ExecuteProcessFinalizeSolutionStep
from KratosMultiphysics.StatisticsApplication.test_utilities import GetInitialVariableValue


class SpatialStatisticsProcessMPITest(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test_model_part")
        self.model_part.SetBufferSize(1)
        self.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)

        input_filename = "spatial_statistics_process/spatial_statistics_process"

        self.__AddNodalSolutionStepVariables()
        communicator = Kratos.DataCommunicator.GetDefault()
        if communicator.Size() == 1:
            self.skipTest("Test can be run only using more than one mpi process")

        if communicator.Rank() == 0 :

            # Original .mdpa file reading
            model_part_io = Kratos.ModelPartIO(input_filename)

            # Partition of the original .mdpa file
            number_of_partitions = communicator.Size() # Number of partitions equals the number of processors
            domain_size = self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
            verbosity = 0
            sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of

            partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
            partitioner.Execute()

            Kratos.Logger.PrintInfo("SpatialStatisticsProcessMPITest","Metis divide finished.")

        communicator.Barrier()

        ## Read the partitioned .mdpa files
        mpi_input_filename = input_filename + "_" + str(communicator.Rank())
        model_part_io = Kratos.ModelPartIO(mpi_input_filename)
        model_part_io.ReadModelPart(self.model_part)

        ## Construct and execute the Parallel fill communicator
        ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()

        InitializeModelPartVariables(self.model_part, False)
        self.model_part.ProcessInfo[Kratos.STEP] = 0

    def tearDown(self):
        rank = self.model_part.GetCommunicator().GetDataCommunicator().Rank()
        if rank == 0:
            kratos_utilities.DeleteFileIfExisting("test_mpi_communicator.time")
        kratos_utilities.DeleteFileIfExisting("spatial_statistics_process/spatial_statistics_process_"+str(rank)+".mdpa")
        kratos_utilities.DeleteFileIfExisting("spatial_statistics_process/spatial_statistics_process_"+str(rank)+".time")
        self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

    def testSpatialStatisticsProcessNodalHistorical(self):
        settings = SpatialStatisticsProcessMPITest.__GetDefaultSettings("nodal_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessNodalNonHistorical(self):
        settings = SpatialStatisticsProcessMPITest.__GetDefaultSettings("nodal_non_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessElementNonHistorical(self):
        settings = SpatialStatisticsProcessMPITest.__GetDefaultSettings("element_non_historical")
        self.__TestMethod(settings)

    def testSpatialStatisticsProcessConditionNonHistorical(self):
        settings = SpatialStatisticsProcessMPITest.__GetDefaultSettings("condition_non_historical")
        self.__TestMethod(settings)


    def __TestMethod(self, settings):
        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        InitializeProcesses(self)

        communicator = self.model_part.GetCommunicator()
        communicator.GetDataCommunicator().Barrier()

        for step in range(0, 12, 2):
            self.model_part.CloneTimeStep(step)
            communicator.GetDataCommunicator().Barrier()
            self.model_part.ProcessInfo[Kratos.STEP] = step
            communicator.GetDataCommunicator().Barrier()
            InitializeModelPartVariables(self.model_part, False)
            ExecuteProcessFinalizeSolutionStep(self)

        for process in self.process_list:
            communicator.GetDataCommunicator().Barrier()
            process.ExecuteFinalize()

    @staticmethod
    def __GetDefaultSettings(container_name):
        # settings_str = r'''
        #     {
        #         "kratos_module" : "KratosMultiphysics.StatisticsApplication",
        #         "python_module" : "spatial_statistics_process",
        #         "Parameters" : {
        #             "model_part_name" : "test_model_part",
        #             "input_variable_settings" : [
        #                 {
        #                     "method_name"    : "median",
        #                     "norm_type"      : "magnitude",
        #                     "container"      : "<CONTAINER_NAME>",
        #                     "variable_names" : ["PRESSURE"],
        #                     "method_settings": {}
        #                 }
        #             ],
        #             "output_settings" : {
        #                 "output_control_variable": "STEP",
        #                 "output_time_interval"   : 1,
        #                 "write_kratos_version"   : false,
        #                 "write_time_stamp"       : false,
        #                 "output_file_settings"   : {
        #                     "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
        #                     "folder_name": "spatial_statistics_process",
        #                     "write_buffer_size" : -1
        #                 }
        #             }
        #         }
        #     }
        # '''
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
        # print(final_settings_str)
        return Kratos.Parameters(final_settings_str)

    def __AddNodalSolutionStepVariables(self):
        # input variables
        self.model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        self.model_part.AddNodalSolutionStepVariable(
            Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        # output variables for output_1
        self.model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_NORM)
        self.model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_NORM)
        self.model_part.AddNodalSolutionStepVariable(Kratos.YIELD_STRESS)
        self.model_part.AddNodalSolutionStepVariable(Kratos.CUTTED_AREA)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.NET_INPUT_MATERIAL)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.WET_VOLUME)


if __name__ == '__main__':
    KratosUnittest.main()