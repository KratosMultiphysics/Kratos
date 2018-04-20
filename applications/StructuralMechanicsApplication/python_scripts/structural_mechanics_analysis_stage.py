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
from base_analysis_stage import AnalysisStage

# Other imports
import sys

class StructuralMechanicsAnalysisStage(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        super(StructuralMechanicsAnalysisStage, self).__init__(model, project_parameters)

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

    def RunSolutionLoop(self):
        while self.time < self.end_time:
            self.time = self.solver.AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

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
        """ Initialize the timestep and advance in time. Called once per timestep """
        super(StructuralMechanicsAnalysisStage, self).InitializeSolutionStep()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("STEP: ", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
            KratosMultiphysics.Logger.PrintInfo("TIME: ", self.time)
        sys.stdout.flush()

        if (self.output_post == True):
            self.gid_output.ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        """ Finalizing the timestep and printing the output. Called once per timestep """
        super(StructuralMechanicsAnalysisStage, self).FinalizeSolutionStep()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalizeSolutionStep()


    def OutputSolutionStep(self):
        """ Printing the output. Called once per timestep """
        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        if (self.output_post == True) and (self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

        self.solver.SaveRestart() # whether a restart-file is written is decided internally

    def Finalize(self):
        super(StructuralMechanicsAnalysisStage, self).Finalize()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalize()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")


    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Structure model part definition
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
        import python_solvers_wrapper_structural
        self.solver = python_solvers_wrapper_structural.CreateSolver(self.main_model_part, self.project_parameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        if not self.using_external_model_part:
            ## Read the model - note that SetBufferSize is done here
            self.solver.ReadModelPart() # TODO move to global instance

    def _InitializeIO(self):
        """ Initialize GiD  I/O """
        self.output_post  = self.project_parameters.Has("output_configuration")
        if (self.output_post == True):
            if (self.parallel_type == "OpenMP"):
                from gid_output_process import GiDOutputProcess as output_process
            elif (self.parallel_type == "MPI"):
                from gid_output_process_mpi import GiDOutputProcessMPI as output_process

            self.gid_output = output_process(self.solver.GetComputingModelPart(),
                                             self.project_parameters["problem_data"]["problem_name"].GetString(),
                                             self.project_parameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

    def _ExecuteInitialize(self):
        """ Initializing the Analysis """
        ## ModelPart is being prepared to be used by the solver
        self.solver.PrepareModelPartForSolver()

        ## Adds the Dofs if they don't exist
        self.solver.AddDofs()

        # Initialize IO
        self._InitializeIO()

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
        self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["loads_process_list"])
        if (self.project_parameters.Has("list_other_processes") == True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["list_other_processes"])
        if (self.project_parameters.Has("json_output_process") == True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_output_process"])
        # Processes for tests
        if (self.project_parameters.Has("json_check_process") == True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_check_process"])
        if (self.project_parameters.Has("check_analytic_results_process") == True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["check_analytic_results_process"])

        if self.is_printing_rank and self.echo_level > 1:
            count = 0
            for process in self.list_of_processes:
                count += 1
                # KratosMultiphysics.Logger.PrintInfo("Process " + str(count), process) # FIXME


    def _ExecuteBeforeSolutionLoop(self):
        """ Perform Operations before the SolutionLoop """
        if (self.output_post == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

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

        if self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True:
            self.time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = start_time
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def GetModelPart(self):
        return self.main_model_part

    def GetSolver(self):
        return self.solver

    def _GetSimulationName(self):
        return "::[KSM Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis_stage.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis_stage.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysisStage(model, parameters)
    simulation.Run()
