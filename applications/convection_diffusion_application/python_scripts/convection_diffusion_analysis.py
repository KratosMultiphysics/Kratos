from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Importing the solvers (if available)
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "not imported")
try:
    import KratosMultiphysics.EigenSolversApplication
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("EigenSolversApplication", "not imported")

# Importing the base class
from analysis_stage import AnalysisStage

# Other imports
import sys

class ConvectionDiffusionAnalysis(AnalysisStage):
    """
    This class is the main-script of the ConvectionDiffusionApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        super(ConvectionDiffusionAnalysis, self).__init__(model, project_parameters)

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            self.is_printing_rank = (KratosMPI.mpi.rank == 0)
        else:
            self.is_printing_rank = True

        self._CreateSolver()

    def Initialize(self):
        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()
        self._ExecuteInitialize()
        self._SetUpListOfProcesses()
        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        ## Solver initialization
        self.solver.Initialize()

        self._ExecuteBeforeSolutionLoop()

    def InitializeSolutionStep(self):
        super(ConvectionDiffusionAnalysis, self).InitializeSolutionStep()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def OutputSolutionStep(self):
        if self.have_output and self.output.IsOutputStep():

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            self.output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

        self.solver.SaveRestart() # whether a restart-file is written is decided internally

    def Finalize(self):
        super(ConvectionDiffusionAnalysis, self).Finalize()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")


    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Convection-diffusion model part definition
        main_model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        if self.model.HasModelPart(main_model_part_name):
            self.main_model_part = self.model[main_model_part_name]
            self.using_external_model_part = True
        else:
            self.main_model_part = KratosMultiphysics.ModelPart(main_model_part_name)
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,
                                                      self.project_parameters["problem_data"]["domain_size"].GetInt())
            self.using_external_model_part = False

        ## Solver construction
        import python_solvers_wrapper_convection_diffusion as solver_wrapper
        self.solver = solver_wrapper.CreateSolver(self.main_model_part, self.project_parameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        if not self.using_external_model_part:
            ## Read the model - note that SetBufferSize is done here
            self.solver.ReadModelPart() # TODO move to global instance

    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        self.have_output = self.project_parameters.Has("output_configuration")
        if self.have_output:
            self.__CheckForDeprecatedGiDSettings()
            if self.parallel_type == "OpenMP":
                from gid_output_process import GiDOutputProcess as OutputProcess
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

            self.output = OutputProcess(self.solver.GetComputingModelPart(),
                                        self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                        self.project_parameters["output_configuration"])

    def _ExecuteInitialize(self):
        """ Initializing the Analysis """
        ## ModelPart is being prepared to be used by the solver
        self.solver.PrepareModelPartForSolver()

        ## Adds the Dofs if they don't exist
        self.solver.AddDofs()

        ## Add the Modelpart to the Model if it is not already there
        if not self.using_external_model_part:
            self.model.AddModelPart(self.main_model_part)

        ## Print model_part and properties
        if self.is_printing_rank and self.echo_level > 1:
            KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
            for properties in self.main_model_part.Properties:
                KratosMultiphysics.Logger.PrintInfo("Property " + str(properties.Id), properties)

    def _SetUpListOfProcesses(self):
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        self.list_of_processes = factory.ConstructListOfProcesses(self.project_parameters["constraints_process_list"])
        self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["fluxes_process_list"])
        if (self.project_parameters.Has("list_other_processes") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["list_other_processes"])
        if (self.project_parameters.Has("json_output_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_output_process"])
        # Processes for tests
        if (self.project_parameters.Has("json_check_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_check_process"])
        if (self.project_parameters.Has("check_analytic_results_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["check_analytic_results_process"])

        #TODO this should be generic
        # initialize GiD  I/O
        self._SetUpGiDOutput()
        if self.have_output:
            self.list_of_processes += [self.output,]

        if self.is_printing_rank and self.echo_level > 1:
            count = 0
            for process in self.list_of_processes:
                count += 1
                # KratosMultiphysics.Logger.PrintInfo("Process " + str(count), process) # FIXME

    def _ExecuteBeforeSolutionLoop(self):
        """ Perform Operations before the SolutionLoop """

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Writing the full ProjectParameters file before solving
        if self.is_printing_rank and self.echo_level > 1:
            f = open("ProjectParametersOutput.json", 'w')
            f.write(self.project_parameters.PrettyPrintJsonString())
            f.close()

        ## Stepping and time settings
        self.solver.SetDeltaTime(self.project_parameters["problem_data"]["time_step"].GetDouble())
        start_time = self.project_parameters["problem_data"]["start_time"].GetDouble()
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] is True:
            self.time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = start_time
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def _GetSimulationName(self):
        return "::[Convection-Diffusion Simulation]:: "

    def __CheckForDeprecatedGiDSettings(self):
        pass

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 convection_diffusion_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 convection_diffusion_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ConvectionDiffusionAnalysis(model, parameters)
    simulation.Run()
