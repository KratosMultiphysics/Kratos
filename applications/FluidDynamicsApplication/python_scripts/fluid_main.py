from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

class FluidMain(object):

    def __init__(self,parameter_file_name='ProjectParameters.json'):

        with open(parameter_file_name,'r') as parameter_file:
            self.project_parameters = Parameters( parameter_file.read() )

        # If this is an MPI run, load the required modules
        
    def SetUpModel(self):
        '''Initialize the model part for the problem (stored as self.model_part) and other general model data.'''

        model_part_name = self.ProjectParameters["problem_data"]["model_part_name"].GetString()
        self.input_model_part = ModelPart(model_part_name)

        self.domain_size = ProjectParameters["problem_data"]["domain_size"].GetInt()
        self.input_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.domain_size)

        #TODO replace this "model" for real one once available
        self.model = { model_part_name : self.input_model_part }

        solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.input_model_part, ProjectParameters["solver_settings"])

        self.solver.AddVariables()
        self.solver.ImportModelPart()

        self.solver.AddDofs()

        self.model_part = self.solver.GetComputingModelPart()

        #TODO: replace MODEL for the Kratos one ASAP
        ##get the list of the skin SubModelParts in the object model
        skin_part_list = self.project_parameters["solver_settings"]["skin_parts"]
        for i in range(skin_part_list.size()):
            part_name = skin_part_list[i].GetString()
            self.model.update( {part_name: main_model_part.GetSubModelPart(part_name)} )


    def SetUpConditions(self):
        '''Read the boundary and initial conditions for the problem and initialize the processes that will manage them.'''

        boundary_condition_process_list = self.project_parameters["boundary_conditions_process_list"]
        self.boundary_condition_processes = list()

        for i in range(boundary_condition_process_list.size()):
            process_settings = boundary_condition_process_list[i]
            module = __import__(process_settings["kratos_module"].GetString())
            interface_file = __import__(process_settings["python_module"].GetString())
            process = interface_file.Factory(process_settings, self.model)
            self.boundary_condition_processes.append( process )

        for process in self.boundary_condition_processes:
            process.ExecuteInitialize()

    def SetUpSolution(self):
        '''Initialize the Python solver and its auxiliary tools and processes.'''

        self.solver.Initialize()

        #TODO this should be generic
        # initialize GiD  I/O
        from gid_output import GiDOutput
        output_settings = ProjectParameters["output_configuration"]
        self.gid_io = GiDOutput(output_settings["output_filename"].GetString(),
                                output_settings["volume_output"].GetBool(),
                                output_settings["gid_post_mode"].GetString(),
                                output_settings["gid_multi_file_flag"].GetString(),
                                output_settings["gid_write_mesh_flag"].GetBool(),
                                output_settings["gid_write_conditions_flag"].GetBool())
        self.output_time = output_settings["output_time"].GetDouble()

        self.gid_io.initialize_results(self.model_part)

        self.nodal_results = []
        for i in range(output_settings["nodal_results"].size()):
            self.nodal_results.append(output_settings["nodal_results"][i].GetString())
        self.gauss_points_results = []
        for i in range(output_settings["gauss_points_results"].size()):
            self.gauss_points_results.append(output_settings["gauss_points_results"][i].GetString())

        for process in self.boundary_condition_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.dt = self.project_parameters["problem_data"]["time_step"].GetDouble()
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

    def Solve(self):
        '''The main solution loop.'''
        
        time = 0.0
        step = 0
        out = 0.0

        while time <= self.end_time:
            time = time + self.dt
            step = step + 1
            
            self.model_part.CloneTimeStep(time)

            print("STEP = ", step)
            print("TIME = ", time)
            
            for process in self.boundary_condition_processes:
                process.ExecuteInitializeSolutionStep()
        
            self.solver.Solve()
        
            # shouldn't this go at the end of the iteration???
            for process in self.boundary_condition_processes:
                process.ExecuteFinalizeSolutionStep()

            if self.output_time <= out:
                for process in self.boundary_condition_processes:
                    process.ExecuteBeforeOutputStep()
    
                self.gid_io.write_results(time,self.model_part,self.nodal_results,self.gauss_points_results)
                out = 0.0
        
                for process in self.boundary_condition_processes:
                    process.ExecuteAfterOutputStep()

            out = out + self.dt

    def FinalizeSolution(self):
        '''Finalize and close open files.'''
        self.gid_io.finalize_results()

        for process in self.boundary_condition_processes:
            process.ExecuteFinalize()

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
