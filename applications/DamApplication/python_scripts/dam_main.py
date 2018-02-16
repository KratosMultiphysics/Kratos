from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
import time as timer
print(timer.ctime())
initial_time = timer.perf_counter()

## Importing modules -----------------------------------------------------------------------------------------

# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
#import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

class Solution(object):

    def LoadParametersFile(self):
        parameter_file = open("ProjectParameters.json",'r')
        self.ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())
        
    def DefineParallelType(self):
        self.parallel_type = self.ProjectParameters["problem_data"]["parallel_type"].GetString()
        parallel=KratosMultiphysics.OpenMPUtils()
        parallel.SetNumThreads(self.ProjectParameters["problem_data"]["number_of_threads"].GetInt())
        if self.parallel_type == "MPI":
            import KratosMultiphysics.mpi as KratosMPI
            print("MPI parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
        else:
            print("OpenMP parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())

    def __init__(self):

        self.LoadParametersFile()
        self.DefineParallelType()
        self.DefineVariables()
        self.CreateModelPart()
        self.SetSolver()

    def DefineVariables(self):
        self.domain_size = ProjectParameters["problem_data"]["domain_size"].GetInt()
        self.problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()
        self.problem_path = os.getcwd()
        self.echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()
        self.buffer_size = ProjectParameters["solver_settings"]["buffer_size"].GetInt()
        self.use_streamline_utility = ProjectParameters["problem_data"]["streamlines_utility"].GetBool()
        self.delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
        self.end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
        self.time = ProjectParameters["problem_data"]["start_time"].GetDouble()
        self.tol = delta_time*1.0e-10
        self.time_scale = ProjectParameters["problem_data"]["time_scale"].GetString()

        # Time Units Converter
        if(self.time_scale=="Weeks"):               # Factor to pass from weeks to seconds
            self.time_unit_converter = 604800.0
        elif(self.time_scale=="Days"):               # Factor to pass from days to seconds
            self.time_unit_converter = 86400.0
        elif(self.time_scale=="Hours"):              # Factor to pass from hours to seconds
            self.time_unit_converter = 3600.0
        else:                                       # No changes
            self.time_unit_converter = 1.0

        # Update time variables
        self.start_time = self.time
        self.delta_time = self.delta_time * self.time_unit_converter
        self.end_time = self.end_time * self.time_unit_converter
        self.time = self.time * self.time_unit_converter
        self.tol = self.tol * self.time_unit_converter

    def CreateModelPart(self):
        self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, time_unit_converter)

    def SetSolver(self):
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])
        
    def Run(self):
        self.Initialize()

        self.RunMainTemporalLoop()

        self.Finalize()

    def Initialize(self):

        self.solver.AddVariables() # Add problem variables
        self.solver.ImportModelPart() # Read model_part (note: the buffer_size is set here)
        self.solver.AddDofs() # Add degrees of freedom

        DamModel = KratosMultiphysics.Model() # Creation of Kratos model
        DamModel.AddModelPart(self.main_model_part)

        # Print model_part and properties
        if(echo_level > 1):
            print(self.main_model_part)
            for self.properties in self.main_model_part.Properties:
                print(self.properties)

        # Construct processes to be applied
        import process_factory
        self.list_of_processes = process_factory.KratosProcessFactory(DamModel).ConstructListOfProcesses( self.ProjectParameters["constraints_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(DamModel).ConstructListOfProcesses( self.ProjectParameters["loads_process_list"] )

        # Print list of constructed processes
        if(echo_level>1):
            for self.process in self.list_of_processes:
                print(self.process)

        # Initialize processes
        for self.process in self.list_of_processes:
            self.process.ExecuteInitialize()

        # Set TIME and DELTA_TIME and fill the previous steps of the buffer with the initial conditions
        self.time = self.time - (buffer_size-1)*self.delta_time
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)
        for step in range(buffer_size-1):
            self.time = self.time + self.delta_time
            self.main_model_part.CloneTimeStep(self.time)

        # Initialize GiD I/O
        computing_model_part = self.solver.GetComputingModelPart()
        output_settings = self.ProjectParameters["output_configuration"]
        if self.parallel_type == "OpenMP":
            import poromechanics_cleaning_utility
            poromechanics_cleaning_utility.CleanPreviousFiles(self.problem_path) # Clean previous post files
            from gid_dam_output_process import GiDDamOutputProcess
            self.gid_output = GiDDamOutputProcess(computing_model_part,
                                             self.problem_name,
                                             self.start_time,
                                             output_settings)
        else:
            from gid_output_process_mpi import GiDOutputProcessMPI
            self.gid_output = GiDOutputProcessMPI(computing_model_part,
                                             self.problem_name,
                                             self.start_time,
                                             output_settings)
        self.gid_output.ExecuteInitialize()

        self.solver.Initialize() # Initialize the solver

        # ExecuteBeforeSolutionLoop
        for self.process in self.list_of_processes:
            self.process.ExecuteBeforeSolutionLoop()

        self.gid_output.ExecuteBeforeSolutionLoop() # Set results when they are written in a single file

        # Initialize streamlines_output_utility
        self.UseStreamlineUtility = False
        if (self.use_streamline_utility == True and self.domain_size==3):
            self.UseStreamlineUtility = True
            import streamlines_output_utility
            self.streamline_utility = streamlines_output_utility.StreamlinesOutputUtility(domain_size)

        if (echo_level > 1):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(self.ProjectParameters.PrettyPrintJsonString())
            f.close()



    def RunMainTemporalLoop(self):

        while( (self.time+self.tol) <= self.end_time ):

            # Update temporal variables
            self.delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            self.time = self.time + self.delta_time
            self.main_model_part.CloneTimeStep(self.time)

            # Update imposed conditions
            for self.process in self.list_of_processes:
                self.process.ExecuteInitializeSolutionStep()

            self.gid_output.ExecuteInitializeSolutionStep()

            self.solver.Solve() # Solve step

            # streamlines_output_utility
            if (self.UseStreamlineUtility== True):
                self.streamline_utility.ComputeOutputStep( main_model_part ,domain_size)

            self.gid_output.ExecuteFinalizeSolutionStep()

            for self.process in self.list_of_processes:
                self.process.ExecuteFinalizeSolutionStep()

            for self.process in self.list_of_processes:
                self.process.ExecuteBeforeOutputStep()

            # Write GiD results
            if self.gid_output.IsOutputStep():
                self.gid_output.PrintOutput()

            for self.process in list_of_processes:
                self.process.ExecuteAfterOutputStep()


    def Initialize(self):
        self.gid_output.ExecuteFinalize() # Finalizing output files

        for self.process in self.list_of_processes:
            self.process.ExecuteFinalize()

        # Finalizing strategy
        if self.parallel_type == "OpenMP":
            self.solver.Clear()

        # Time control
        print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - initial_time)," seconds.")
        print(timer.ctime())
        
if __name__ == "__main__":
    Solution().Run()
