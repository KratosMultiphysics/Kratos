import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MeshMovingApplication as KMM
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
import os

class MeshMovingAnalysisForTesting(MeshMovingAnalysis):

    def __init__(self, model, project_parameters, mesh_vel_calc_helper):
        super(MeshMovingAnalysisForTesting, self).__init__(model, project_parameters)
        self.mesh_vel_calc_helper = mesh_vel_calc_helper
        # Variables needed for the mesh-velocity-calculation
        self.model["MainModelPart"].AddNodalSolutionStepVariable(KM.MESH_VELOCITY)
        self.model["MainModelPart"].AddNodalSolutionStepVariable(KM.MESH_ACCELERATION)

    def RunSolutionLoop(self):
        """custom to also compute the mesh-velocities"""
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()

            KMM.CalculateMeshVelocities(
                self._GetSolver().GetComputingModelPart(),
                self.mesh_vel_calc_helper)

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

class MeshMovingTestCase(KratosUnittest.TestCase):

    def executeTest(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open('generic_rectangle_test_parameters.json', 'r') as parameter_file:
                self.project_parameters = KM.Parameters(parameter_file.read())

            self.__SetProblemData()
            self.__SetSolverSettings()
            self.__SetPrintOfReferenceResults()
            self.__SetOutputConfiguration()

            self.__SetLoggerSeverity()

            model = KM.Model()
            MeshMovingAnalysisForTesting(model, self.project_parameters, self.mesh_vel_calc_helper).Run()

    def __SetProblemData(self):
        problem_data = self.project_parameters["problem_data"]

        problem_data["problem_name"].SetString(self.__GetProblemName())

        if KM.DataCommunicator.GetDefault().IsDistributed(): # whether this is an mpi-execution
            problem_data["parallel_type"].SetString("MPI")
        else:
            problem_data["parallel_type"].SetString("OpenMP")

    def __SetSolverSettings(self):
        solver_settings = self.project_parameters["solver_settings"]

        solver_settings["domain_size"].SetInt(self.domain_size)
        solver_settings["buffer_size"].SetInt(KM.GetMinimumBufferSize(self.mesh_vel_calc_helper))
        solver_settings["solver_type"].SetString(self.solver_type)
        input_file_name = os.path.join("test_mdpa_files", "rectangle_{}_test".format(self.__GetElementTopology()))
        solver_settings["model_import_settings"]["input_filename"].SetString(input_file_name)

    def __SetPrintOfReferenceResults(self):
        processes = self.project_parameters["processes"]
        result_file_name = os.path.join("test_results", self.__GetProblemName()+"_results.json")
        if self.print_reference_results:
            warn_msg  = 'The test "{}" is configured for writing reference results\n'.format(self.__GetProblemName())
            warn_msg += 'This is only for setting up the test, the results are NOT checked!'
            KM.Logger.PrintWarning("WARNING MeshMovingTestCase", warn_msg)
            processes.AddValue("json_output", KM.Parameters("""[{
                "python_module" : "json_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "JsonOutputProcess",
                "Parameters"    : {
                    "output_variables" : ["MESH_DISPLACEMENT_X",
                                          "MESH_DISPLACEMENT_Y",
                                          "MESH_DISPLACEMENT_Z",
                                          "MESH_VELOCITY_X",
                                          "MESH_VELOCITY_Y",
                                          "MESH_VELOCITY_Z"],
                    "output_file_name" : \""""+result_file_name+"""\",
                    "model_part_name"  : "Probe_1",
                    "time_frequency"   : 0.1,
                    "use_node_coordinates" : true
                }
            }]"""))
        else:
            processes.AddValue("json_check_process", KM.Parameters("""[{
                "python_module" : "from_json_check_result_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "FromJsonCheckResultProcess",
                "Parameters"    : {
                    "check_variables" : ["MESH_DISPLACEMENT_X",
                                         "MESH_DISPLACEMENT_Y",
                                         "MESH_DISPLACEMENT_Z",
                                         "MESH_VELOCITY_X",
                                         "MESH_VELOCITY_Y",
                                         "MESH_VELOCITY_Z"],
                    "input_file_name"  : \""""+result_file_name+"""\",
                    "model_part_name"  : "Probe_1",
                    "time_frequency"   : 0.1,
                    "use_node_coordinates" : true
                }
            }]"""))

    def __SetOutputConfiguration(self):
        output_processes = self.project_parameters["output_processes"]
        problem_name = self.__GetProblemName()
        # to check the results: add output settings block if needed
        if self.print_vtk_output:
            output_processes.AddValue("vtk_output", KM.Parameters("""[{
                "python_module" : "vtk_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "VtkOutputProcess",
                "help"          : "This process writes postprocessing files for Paraview",
                "Parameters"    : {
                    "model_part_name"                    : "MainModelPart",
                    "output_control_type"                : "step",
                    "output_frequency"                   : 1,
                    "file_format"                        : "binary",
                    "output_precision"                   : 7,
                    "output_sub_model_parts"             : false,
                    "folder_name"                        : \""""+problem_name+"""\",
                    "save_output_files_in_folder"        : true,
                    "nodal_solution_step_data_variables" : ["MESH_DISPLACEMENT","MESH_VELOCITY"],
                    "nodal_data_value_variables"         : [],
                    "element_data_value_variables"       : [],
                    "condition_data_value_variables"     : []
                }
            }]"""))

        if self.print_gid_output:
            output_processes.AddValue("gid_output", KM.Parameters("""[{
                "python_module" : "gid_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "GiDOutputProcess",
                "help"          : "This process writes postprocessing files for GiD",
                "Parameters"    : {
                    "model_part_name"        : "MainModelPart",
                    "output_name"            : \""""+problem_name+"""\",
                    "postprocess_parameters" : {
                        "result_file_configuration" : {
                            "gidpost_flags"       : {
                                "GiDPostMode"           : "GiD_PostBinary",
                                "WriteDeformedMeshFlag" : "WriteDeformed",
                                "WriteConditionsFlag"   : "WriteConditions",
                                "MultiFileFlag"         : "SingleFile"
                            },
                            "file_label"          : "time",
                            "output_control_type" : "step",
                            "output_frequency"    : 1.0,
                            "body_output"         : true,
                            "nodal_results"       : ["MESH_DISPLACEMENT","MESH_VELOCITY"]
                        }
                    }
                }
            }]"""))

    def __GetElementTopology(self):
        return "{}D{}N".format(self.domain_size, self.number_of_nodes_per_elements)

    def __GetProblemName(self):
        return "mesh_moving_rectangle_" + self.solver_type + "_" + self.__GetElementTopology()

    def __SetLoggerSeverity(self):
        if self.print_logger_info:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
