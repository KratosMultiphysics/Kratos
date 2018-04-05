from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShallowWaterApplication as Shallow

# Importing the solvers (if available)
try:
    import Kratos.ExternalSolversApplication
    Kratos.Logger.PrintInfo("ExternalSolversApplication", "successfully imported")
except ImportError:
    Kratos.Logger.PrintInfo("ExternalSolversApplication", "not imported")

class ShallowWaterAnalysis(object): # TODO in the future this could derive from a BaseClass in the Core
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, project_parameters):
        if (type(project_parameters) == str): # a file name is provided
            with open(project_parameters,'r') as parameter_file:
                self.ProjectParameters = Kratos.Parameters(parameter_file.read())
        elif (type(project_parameters) == Kratos.Parameters): # a Parameters object is provided
            self.ProjectParameters = project_parameters
        else:
            raise Exception("Input is expected to be provided as a Kratos Parameters object or a file name")
        self.is_printing_rank = True

    #### Public functions to run the Analysis ####
    def Run(self):
        '''Wrapper function for the solution.'''
        self.InitializeAnalysis()
        self.RunMainTemporalLoop()
        self.FinalizeAnalysis()

    def InitializeAnalysis(self):
        '''Wrapper function comprising the definition of the model and the initialization of the problem.'''
        self.SetUpModel()
        self.SetUpProcesses()
        self.SetUpAnalysis()

    def RunMainTemporalLoop(self):
        '''The main solution loop.'''
        while self.time <= self.end_time:
            dt = self.solver.ComputeDeltaTime()
            self.time = self.time + dt
            self.step = self.step + 1

            self.main_model_part.CloneTimeStep(self.time)
            self.main_model_part.ProcessInfo[Kratos.STEP] = self.step

            if self.is_printing_rank:
                Kratos.Logger.PrintInfo("Fluid Dynamics Analysis","STEP = ", self.step)
                Kratos.Logger.PrintInfo("Fluid Dynamics Analysis","TIME = ", self.time)

            self.InitializeSolutionStep()
            self.SolveSingleStep()
            self.FinalizeSolutionStep()

    def FinalizeAnalysis(self):
        '''Finalize the simulation and close open files.'''
        for process in self.list_of_processes:
            process.ExecuteFinalize()
        if self.have_output:
            self.output.ExecuteFinalize()

    def SetUpModel(self):
        '''Initialize the model part for the problem and other general model data.'''
        
        ## Defining variables ----------------------------------------------------------------------------------------
        self.problem_name   = self.ProjectParameters["problem_data"]["problem_name"].GetString()
        self.parallel_type  = self.ProjectParameters["problem_data"]["parallel_type"].GetString()
        self.echo_level     = self.ProjectParameters["solver_settings"]["echo_level"].GetInt()
        self.end_time       = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
        self.time           = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        self.step           = 0
        domain_size         = self.ProjectParameters["problem_data"]["domain_size"].GetInt()
        gravity             = self.ProjectParameters["problem_data"]["gravity"].GetDouble()
        time_scale          = self.ProjectParameters["problem_data"]["time_scale"].GetString()
        water_height_scale  = self.ProjectParameters["problem_data"]["water_height_scale"].GetString()

        # Time unit converter
        if   time_scale == "seconds":
            time_unit_converter =     1
        elif time_scale == "minutes":
            time_unit_converter =    60
        elif time_scale == "hours":
            time_unit_converter =  3600
        elif time_scale == "days":
            time_unit_converter = 86400
        else:
            raise Exception("unknown time scale")

        # Water height unit converter
        if   water_height_scale == "meters":
            water_height_unit_converter = 1.0
        elif water_height_scale == "millimeters":
            water_height_unit_converter = 0.001
        else:
            raise Exception("unknown water height scale")

        ## Model part ------------------------------------------------------------------------------------------------

        # Defining the model part
        self.main_model_part = Kratos.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, domain_size)
        self.main_model_part.ProcessInfo.SetValue(Kratos.TIME, self.time)
        self.main_model_part.ProcessInfo.SetValue(Kratos.GRAVITY_Z, gravity * time_unit_converter**2)
        self.main_model_part.ProcessInfo.SetValue(Shallow.TIME_UNIT_CONVERTER, time_unit_converter)
        self.main_model_part.ProcessInfo.SetValue(Shallow.WATER_HEIGHT_UNIT_CONVERTER, water_height_unit_converter)

        # Solver construction (main settings methods are located in the solver_module)
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

        self.solver.AddVariables()
        self.solver.ImportModelPart()
        self.solver.AddDofs()

        # Fill a Model instance using input
        self.model = Kratos.Model()
        self.model.AddModelPart(self.main_model_part)

    def SetUpProcesses(self):
        '''
        Read the definition of initial and boundary conditions for the problem and initialize the processes that will manage them.
        Also initialize any additional processes present in the problem (such as those used to calculate additional results).
        '''
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        # "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: bathymetry is firstly constructed. Initial conditions might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        self.list_of_processes  = factory.ConstructListOfProcesses( self.ProjectParameters["bathymetry_process_list"] )
        self.list_of_processes += factory.ConstructListOfProcesses( self.ProjectParameters["initial_conditions_process_list"] )
        self.list_of_processes += factory.ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )

    def SetUpAnalysis(self):
        '''
        Initialize the Python solver and its auxiliary tools and processes.
        This function should prepare everything so that the simulation
        can start immediately after exiting it.
        '''

        for process in self.list_of_processes:
            process.ExecuteInitialize()

        self.solver.Initialize()

        self._SetUpGiDOutput()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        if self.have_output:
            self.output.ExecuteBeforeSolutionLoop()

    def InitializeSolutionStep(self):

        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        if self.have_output:
            self.output.ExecuteInitializeSolutionStep()

    def SolveSingleStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.Solve()

    def FinalizeSolutionStep(self):

        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        if self.have_output:
            self.output.ExecuteFinalizeSolutionStep()

        if self.have_output and self.output.IsOutputStep():

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            self.output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.step + 1 >= self.solver.GetMinimumBufferSize()

    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        self.have_output = self.ProjectParameters.Has("output_configuration")
        if self.have_output:
            if self.parallel_type == "OpenMP":
                from gid_output_process import GiDOutputProcess as OutputProcess
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

            self.output = OutputProcess(self.main_model_part,
                                        self.ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                        self.ProjectParameters["output_configuration"])

            self.output.ExecuteInitialize()



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

    ShallowWaterAnalysis(project_parameters_file_name).Run()