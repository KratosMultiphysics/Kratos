from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteEmbeddedTest(KratosUnittest.TestCase):

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        self.main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

        Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        ## Solver construction
        import python_solvers_wrapper_fluid
        self.solver = python_solvers_wrapper_fluid.CreateSolver(self.main_model_part, ProjectParameters)

        self.solver.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        self.solver.ImportModelPart()

        ## Add AddDofs
        self.solver.AddDofs()

        ## Initialize GiD  I/O
        self.output_flag = False
        if (self.output_flag == True):
            from gid_output_process import GiDOutputProcess
            self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
                                               ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                               ProjectParameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

        ## Solver initialization
        self.solver.Initialize()

        ## Get the list of the skin submodel parts in the object Model
        for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
            Model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

        ## Get the gravity submodel part in the object Model
        for i in range(ProjectParameters["gravity"].size()):
            gravity_part_name = ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
            Model.update({gravity_part_name: self.main_model_part.GetSubModelPart(gravity_part_name)})

        ## Processes construction
        import process_factory
        self.list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["gravity"] )
        self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()


    def Solve(self):

        ## Stepping and time settings
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = 0.0
        step = 0

        if (self.output_flag == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            Dt = self.solver.ComputeDeltaTime()
            time = time + Dt
            step = step + 1
            self.main_model_part.CloneTimeStep(time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()

            if (self.output_flag == True):
                self.gid_output.ExecuteInitializeSolutionStep()

            if(step >= 3):
                self.solver.Solve()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

            if (self.output_flag == True):
                self.gid_output.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            if (self.output_flag == True):
                self.gid_output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
