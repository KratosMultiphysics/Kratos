from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the solvers (if available)
try:
    import KratosMultiphysics.ExternalSolversApplication
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "not imported")
try:
    import KratosMultiphysics.EigenSolversApplication
    KratosMultiphysics.Logger.PrintInfo("EigenSolversApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("EigenSolversApplication", "not imported")

# Importing the base class
from analysis_stage import AnalysisStage

# Other imports
import sys

class StructuralMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        super(StructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def OutputSolutionStep(self):
        super(StructuralMechanicsAnalysis, self).OutputSolutionStep()

        self._GetSolver().SaveRestart()


    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        import python_solvers_wrapper_structural
        return python_solvers_wrapper_structural.CreateSolver(self.main_model_part, self.project_parameters)

    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance'''
        self.__CheckForDeprecatedGiDSettings()
        if self.parallel_type == "OpenMP":
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        gid_output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                   self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                   self.project_parameters["output_configuration"])

        return gid_output


    # def _ExecuteInitialize(self):
    #     """ Initializing the Analysis """
    #     ## ModelPart is being prepared to be used by the solver
    #     self.solver.PrepareModelPartForSolver()

    #     ## Adds the Dofs if they don't exist
    #     self.solver.AddDofs()

    #     ## Add the Modelpart to the Model if it is not already there
    #     if not self.using_external_model_part:
    #         self.model.AddModelPart(self.main_model_part)

    #     ## Print model_part and properties
    #     if self.is_printing_rank and self.echo_level > 1:
    #         KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
    #         for properties in self.main_model_part.Properties:
    #             KratosMultiphysics.Logger.PrintInfo("Property " + str(properties.Id), properties)

    # def _SetUpListOfProcesses(self):
    #     from process_factory import KratosProcessFactory
    #     factory = KratosProcessFactory(self.model)
    #     self.list_of_processes = factory.ConstructListOfProcesses(self.project_parameters["constraints_process_list"])
    #     self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["loads_process_list"])
    #     if (self.project_parameters.Has("list_other_processes") is True):
    #         self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["list_other_processes"])
    #     if (self.project_parameters.Has("json_output_process") is True):
    #         self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_output_process"])
    #     # Processes for tests
    #     if (self.project_parameters.Has("json_check_process") is True):
    #         self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_check_process"])
    #     if (self.project_parameters.Has("check_analytic_results_process") is True):
    #         self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["check_analytic_results_process"])

    #     #TODO this should be generic
    #     # initialize GiD  I/O
    #     self._SetUpGiDOutput()
    #     if self.have_output:
    #         self.list_of_processes += [self.output,]

    #     if self.is_printing_rank and self.echo_level > 1:
    #         count = 0
    #         for process in self.list_of_processes:
    #             count += 1
    #             # KratosMultiphysics.Logger.PrintInfo("Process " + str(count), process) # FIXME

    # def _ExecuteBeforeSolutionLoop(self):
    #     """ Perform Operations before the SolutionLoop """

    #     for process in self.list_of_processes:
    #         process.ExecuteBeforeSolutionLoop()

    #     ## Writing the full ProjectParameters file before solving
    #     if self.is_printing_rank and self.echo_level > 1:
    #         f = open("ProjectParametersOutput.json", 'w')
    #         f.write(self.project_parameters.PrettyPrintJsonString())
    #         f.close()

    #     ## Stepping and time settings
    #     self.solver.SetDeltaTime(self.project_parameters["problem_data"]["time_step"].GetDouble())
    #     start_time = self.project_parameters["problem_data"]["start_time"].GetDouble()
    #     self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

    #     if self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] is True:
    #         self.time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
    #     else:
    #         self.time = start_time
    #         self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

    #     if self.is_printing_rank:
    #         KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is temporary to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(StructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        if len(list_of_processes) == 0:
            KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsAnalysis", "Using the old way to create the processes, this will be removed!")
            from process_factory import KratosProcessFactory
            factory = KratosProcessFactory(self.model)
            list_of_processes = factory.ConstructListOfProcesses(self.project_parameters["constraints_process_list"])
            list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["loads_process_list"])
            if (self.project_parameters.Has("list_other_processes") is True):
                list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["list_other_processes"])
            if (self.project_parameters.Has("json_output_process") is True):
                list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_output_process"])
            # Processes for tests
            if (self.project_parameters.Has("json_check_process") is True):
                list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_check_process"])
            if (self.project_parameters.Has("check_analytic_results_process") is True):
                list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["check_analytic_results_process"])
        else:
            if self.project_parameters.Has("constraints_process_list"):
                raise Exception("Mixing of process initialization is not alowed!")
            if self.project_parameters.Has("loads_process_list"):
                raise Exception("Mixing of process initialization is not alowed!")
            if self.project_parameters.Has("list_other_processes"):
                raise Exception("Mixing of process initialization is not alowed!")
            if self.project_parameters.Has("json_output_process"):
                raise Exception("Mixing of process initialization is not alowed!")
            if self.project_parameters.Has("json_check_process"):
                raise Exception("Mixing of process initialization is not alowed!")
            if self.project_parameters.Has("check_analytic_results_process"):
                raise Exception("Mixing of process initialization is not alowed!")

        return list_of_processes

    def _GetSimulationName(self):
        return "::[KSM Simulation]:: "

    def __CheckForDeprecatedGiDSettings(self):
        if self.project_parameters["output_configuration"].Has("result_file_configuration"):
            res_file_config = self.project_parameters["output_configuration"]["result_file_configuration"]
            if res_file_config.Has("nodal_results"):
                nodal_res = res_file_config["nodal_results"]
                for i in range(nodal_res.size()):
                    var_name = nodal_res[i].GetString()
                    if var_name == "TORQUE":
                        err_msg  = 'Requesting output for "TORQUE" which is not available any more\n'
                        err_msg += 'It was renamed to "MOMENT_REACTION"'
                        raise Exception(err_msg)

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
