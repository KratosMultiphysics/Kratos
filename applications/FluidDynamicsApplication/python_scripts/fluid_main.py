from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
try:
    from KratosMultiphysics.ExternalSolversApplication import *
except ImportError:
    pass

class FluidMain(object):

    def __init__(self,parameter_file_name='ProjectParameters.json'):

        with open(parameter_file_name,'r') as parameter_file:
            self.project_parameters = Parameters( parameter_file.read() )

        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        # If this is an MPI run, load the distributed memory modules
        if (self.parallel_type == "MPI"):
            from KratosMultiphysics.mpi import mpi
            import KratosMultiphysics.MetisApplication as Metis
            import KratosMultiphysics.TrilinosApplication as Trilinos
            self.is_printing_rank = (mpi.rank == 0)
        else:
            self.is_printing_rank = True

        
    def SetUpModel(self):
        '''Initialize the model part for the problem (stored as self.model_part) and other general model data.'''

        model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        self.main_model_part = ModelPart(model_part_name)

        self.domain_size = self.project_parameters["problem_data"]["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.domain_size)

        ## Solver construction
        import python_solvers_wrapper_fluid
        self.solver = python_solvers_wrapper_fluid.CreateSolver(self.main_model_part, self.project_parameters)

        self.solver.AddVariables()
        self.solver.ImportModelPart()
        self.solver.AddDofs()

        self.model_part = self.solver.GetComputingModelPart()

        # Fill a Model instance using input
        self.model = Model()
        self.model.AddModelPart(self.main_model_part)

        # Add the skin SubModelParts to the model
        for i in range(self.project_parameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = self.project_parameters["solver_settings"]["skin_parts"][i].GetString()
            self.model.AddModelPart(self.main_model_part.GetSubModelPart(skin_part_name))

        # Add the no-skin SubModelParts parts to the model (results processes and no-skin conditions)
        for i in range(self.project_parameters["solver_settings"]["no_skin_parts"].size()):
            no_skin_part_name = self.project_parameters["solver_settings"]["no_skin_parts"][i].GetString()
            self.model.AddModelPart(self.model_part.GetSubModelPart(no_skin_part_name))
            
        # Add the initial conditions SubModelParts to the model
        for i in range(self.project_parameters["initial_conditions_process_list"].size()):
            initial_cond_part_name = self.project_parameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
            self.model.AddModelPart(self.main_model_part.GetSubModelPart(initial_cond_part_name))

        # Add the gravity SubModelParts to the model
        for i in range(self.project_parameters["gravity"].size()):
            gravity_part_name = self.project_parameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
            self.model.AddModelPart(self.main_model_part.GetSubModelPart(gravity_part_name))

    def SetUpConditions(self):
        '''Read the boundary and initial conditions for the problem and initialize the processes that will manage them.'''

        ## Processes construction
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        self.simulation_processes =  factory.ConstructListOfProcesses( self.project_parameters["gravity"] )
        self.simulation_processes += factory.ConstructListOfProcesses( self.project_parameters["initial_conditions_process_list"] )
        self.simulation_processes += factory.ConstructListOfProcesses( self.project_parameters["boundary_conditions_process_list"] )
        self.simulation_processes += factory.ConstructListOfProcesses( self.project_parameters["auxiliar_process_list"] )

        for process in self.simulation_processes:
            process.ExecuteInitialize()

    def SetUpSolution(self):
        '''Initialize the Python solver and its auxiliary tools and processes.'''

        self.solver.Initialize()

        #TODO this should be generic
        # initialize GiD  I/O
        self._SetUpGiDOutput()

        for process in self.simulation_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        ## Writing the full ProjectParameters file before solving
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())


    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        self.have_output = self.project_parameters.Has("output_configuration")
        if self.have_output:
            if self.parallel_type == "OpenMP":
                from gid_output_process import GiDOutputProcess
                self.output = GiDOutputProcess(self.model_part,
                                               self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                               self.project_parameters["output_configuration"])
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI
                self.output = GiDOutputProcessMPI(self.model_part,
                                                  self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                                  self.project_parameters["output_configuration"])

            self.output.ExecuteInitialize()

    def Solve(self):
        '''The main solution loop.'''
        
        time = 0.0
        step = 0

        while time <= self.end_time:
            dt = self.solver.ComputeDeltaTime()
            time = time + dt
            step = step + 1
            
            self.main_model_part.CloneTimeStep(time)
            self.main_model_part.ProcessInfo[STEP] = step

            if self.is_printing_rank:
                Logger.Print("STEP = ", step)
                Logger.Print("TIME = ", time)
            
            for process in self.simulation_processes:
                process.ExecuteInitializeSolutionStep()

            if self.have_output:
                self.output.ExecuteInitializeSolutionStep()
        
            self.solver.Solve()
        
            # shouldn't this go at the end of the iteration???
            for process in self.simulation_processes:
                process.ExecuteFinalizeSolutionStep()

            if self.have_output:
                self.output.ExecuteFinalizeSolutionStep()

            if self.have_output and self.output.IsOutputStep():

                for process in self.simulation_processes:
                    process.ExecuteBeforeOutputStep()
    
                self.output.PrintOutput()
        
                for process in self.simulation_processes:
                    process.ExecuteAfterOutputStep()

    def FinalizeSolution(self):
        '''Finalize and close open files.'''

        for process in self.simulation_processes:
            process.ExecuteFinalize()

        if self.have_output:
            self.output.ExecuteFinalize()

    def Run(self):
        '''Wrapper function for the solution.'''
        self.SetUpModel()
        self.SetUpConditions()
        self.SetUpSolution()
        self.Solve()
        self.FinalizeSolution()

if __name__ == '__main__':
    solver = FluidMain()
    solver.Run()
