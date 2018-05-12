from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

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


class MeshMovingAnalysis(object): # TODO in the future this could derive from a BaseClass in the Core
    """
    This class is the main-script of the MeshMovingApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, ProjectParameters, external_model_part=None):
        if (type(ProjectParameters) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")
        self.ProjectParameters = ProjectParameters
        self.__CreateSolver(external_model_part)

    #### Public functions to run the Analysis ####
    def Run(self):
        # Attention, it does not make too much sense calling this function, since MeshMotion is usually coupled with another solver
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

    def RunMainTemporalLoop(self):
        # Attention, it does not make too much sense calling this function, since MeshMotion is usually coupled with another solver
        while self.time < self.end_time:
            self.InitializeTimeStep()
            self.SolveTimeStep()
            self.FinalizeTimeStep()

    #### Public functions defining the Interface to the CoSimulationApplication ####
    def Initialize(self):
        """ This function must be called only once! """
        self.__ExecuteInitialize()
        self.__InitializeIO()
        self.__ExecuteBeforeSolutionLoop()

    def InitializeTimeStep(self):
        """ This function must be called once at the beginning of each timestep """
        self.__ExecuteInitializeSolutionStep()

    def SolveTimeStep(self):
        """ This function can be called several times within each timestep """
        self.__SolveSolutionStep()

    def FinalizeTimeStep(self):
        """ This function must be called once at the end of each timestep """
        self.__ExecuteFinalizeSolutionStep()

    def Finalize(self):
        """ This function must be called only once! """
        self.__ExecuteFinalize()

    #### Internal functions ####
    def __CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not passed from outside) """
        if external_model_part != None:
            # This is a temporary solution until the importing of the ModelPart
            # is removed from the solver (needed e.g. for Optimization)
            if (type(external_model_part) != KratosMultiphysics.ModelPart):
                raise Exception("Input is expected to be provided as a Kratos ModelPart object")
            self.using_external_model_part = True
        else:
            self.using_external_model_part = False

        ## Get echo level and parallel type
        self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.ProjectParameters["problem_data"]["parallel_type"].GetString()

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            self.is_printing_rank = (KratosMPI.mpi.rank == 0)
        else:
            self.is_printing_rank = True

        ## Structure model part definition
        if self.using_external_model_part:
            self.main_model_part = external_model_part
        else:
            main_model_part_name = self.ProjectParameters["problem_data"]["model_part_name"].GetString()
            self.main_model_part = KratosMultiphysics.ModelPart(main_model_part_name)
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,
                                                      self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        ## Solver construction
        import python_solvers_wrapper_mesh_motion
        self.solver = python_solvers_wrapper_mesh_motion.CreateSolver(self.main_model_part, self.ProjectParameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        if not self.using_external_model_part:
            ## Read the model - note that SetBufferSize is done here
            self.solver.ReadModelPart() # TODO move to global instance

    def __InitializeIO(self):
        """ Initialize GiD  I/O """
        self.output_post  = self.ProjectParameters.Has("output_configuration")
        if (self.output_post == True):
            if (self.parallel_type == "OpenMP"):
                from gid_output_process import GiDOutputProcess as output_process
            elif (self.parallel_type == "MPI"):
                from gid_output_process_mpi import GiDOutputProcessMPI as output_process

            self.gid_output = output_process(self.solver.GetComputingModelPart(),
                                             self.ProjectParameters["problem_data"]["problem_name"].GetString(),
                                             self.ProjectParameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

    def __ExecuteInitialize(self):
        """ Initializing the Analysis """
        ## ModelPart is being prepared to be used by the solver
        self.solver.PrepareModelPartForSolver()

        ## Adds the Dofs if they don't exist
        self.solver.AddDofs()

        ## Creation of the Kratos model (build sub_model_parts or submeshes)
        self.model = KratosMultiphysics.Model()
        self.model.AddModelPart(self.main_model_part)

        ## Print model_part and properties
        if self.is_printing_rank and (self.echo_level > 1):
            KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)

        ## Processes construction
        import process_factory
        self.list_of_processes = []
        if (self.ProjectParameters.Has("boundary_conditions_process_list") == True):
            self.list_of_processes = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["boundary_conditions_process_list"])
        if (self.ProjectParameters.Has("list_other_processes") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])
        if (self.ProjectParameters.Has("json_output_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["json_output_process"])
        # Processes for tests
        if (self.ProjectParameters.Has("json_check_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["json_check_process"])
        if (self.ProjectParameters.Has("check_analytic_results_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["check_analytic_results_process"])

        if self.is_printing_rank and (self.echo_level > 1):
            count = 0
            for process in self.list_of_processes:
                count += 1
                # KratosMultiphysics.Logger.PrintInfo("Process " + str(count), process) # FIXME

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        ## Solver initialization
        self.solver.Initialize()

    def __ExecuteBeforeSolutionLoop(self):
        """ Perform Operations before the SolutionLoop """
        if (self.output_post == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Writing the full ProjectParameters file before solving
        if self.is_printing_rank and (self.echo_level > 1):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(self.ProjectParameters.PrettyPrintJsonString())
            f.close()

        ## Stepping and time settings
        self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        self.end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        if self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True:
            self.time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = start_time
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -START- ")

    def __ExecuteInitializeSolutionStep(self):
        """ Initialize the timestep and advance in time. Called once per timestep """
        self.time += self.delta_time
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(self.time)

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("STEP: ", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
            KratosMultiphysics.Logger.PrintInfo("TIME: ", self.time)

        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        if (self.output_post == True):
            self.gid_output.ExecuteInitializeSolutionStep()

    def __ExecuteBeforeSolve(self):
        """ Function to be called before solving. Can be executed several times per timestep """
        pass

    def __SolveSolutionStep(self):
        """ Solving one step. Can be called several times per timestep """
        self.__ExecuteBeforeSolve()
        self.solver.Solve()
        self.__ExecuteAfterSolve()

    def __ExecuteAfterSolve(self):
        """ Function to be called after solving. Can be executed several times per timestep """
        pass

    def __ExecuteFinalizeSolutionStep(self):
        """ Finalizing the timestep and printing the output. Called once per timestep """
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalizeSolutionStep()

        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        if (self.output_post == True) and (self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

    def __ExecuteFinalize(self):
        """ Operations to be performed at the end of the Analysis """
        for process in self.list_of_processes:
            process.ExecuteFinalize()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalize()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -END- ")

    def GetModelPart(self):
        return self.main_model_part

    def GetSolver(self):
        return self.solver


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 mesh_moving_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 mesh_moving_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

    MeshMovingAnalysis(ProjectParameters).Run()