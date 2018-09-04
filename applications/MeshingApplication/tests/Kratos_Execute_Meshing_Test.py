from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
try:
  import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
  missing_external_fluid_dependencies = False
except ImportError as e:
    missing_external_fluid_dependencies = True
try:
  import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
  missing_external_solid_dependencies = False
except ImportError as e:
    missing_external_solid_dependencies = True

import os
import process_factory

class Kratos_Execute_Test:

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        # To avoid many prints
        echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()
        if (echo_level == 0):
            KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)


        self.model = KratosMultiphysics.Model()
        self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        self.problem_type = self.ProjectParameters["problem_data"]["problem_type"].GetString()
        self.solve_problem = self.ProjectParameters["problem_data"]["solve_problem"].GetBool()

        if (self.problem_type  == "fluid" and missing_external_fluid_dependencies == False):
            ## Solver construction
            import python_solvers_wrapper_fluid as fluid_wrapper
            self.solver = fluid_wrapper.CreateSolver(self.model, self.ProjectParameters)
        elif (self.problem_type  == "solid" and missing_external_solid_dependencies == False):
            self.ProjectParameters["solver_settings"].AddValue("model_part_name", self.ProjectParameters["problem_data"]["model_part_name"])
            self.ProjectParameters["solver_settings"].AddValue("domain_size", self.ProjectParameters["problem_data"]["domain_size"])
            time_stepping_params = KratosMultiphysics.Parameters("{}")
            time_stepping_params.AddValue("time_step", self.ProjectParameters["problem_data"]["time_step"])
            self.ProjectParameters["solver_settings"].AddValue("time_stepping", time_stepping_params)
            params = KratosMultiphysics.Parameters("""{"computing_model_part_name": "Structure"}""")
            self.ProjectParameters["solver_settings"].AddValue("computing_model_part_name", params["computing_model_part_name"])
            # Construct the solver (main setting methods are located in the solver_module)
            import python_solvers_wrapper_structural
            self.solver = python_solvers_wrapper_structural.CreateSolver(self.model, self.ProjectParameters)
        else:
            raise NameError('Problem type not defined or failing in the import')

        # Add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
        # If we integrate it in the model part we cannot use combined solvers
        self.solver.AddVariables()

        # Read model_part (note: the buffer_size is set here) (restart can be read here)
        self.solver.ImportModelPart()
        self.solver.PrepareModelPart()

        # Add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
        # If we integrate it in the model part we cannot use combined solvers
        self.solver.AddDofs()

        # ### Output settings start ####
        self.problem_path = os.getcwd()
        self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()

        # ### Output settings start ####
        self.output_post = ProjectParameters.Has("output_configuration")
        if (self.output_post == True):
            from gid_output_process import GiDOutputProcess
            output_settings = ProjectParameters["output_configuration"]
            self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
                                               self.problem_name,
                                               output_settings)
            self.gid_output.ExecuteInitialize()

        # Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        self.solver.SetEchoLevel(0) # Avoid to print anything

        if self.problem_type  == "fluid":
            # Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
            # #Get the list of the submodel part in the object Model
            for i in range(self.ProjectParameters["solver_settings"]["skin_parts"].size()):
                skin_part_name = self.ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
                self.model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

            ## Get the list of the initial conditions submodel parts in the object Model
            for i in range(self.ProjectParameters["initial_conditions_process_list"].size()):
                initial_cond_part_name = self.ProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
                self.model.update({initial_cond_part_name: self.main_model_part.GetSubModelPart(initial_cond_part_name)})

            ## Get the gravity submodel part in the object Model
            for i in range(self.ProjectParameters["gravity"].size()):
                gravity_part_name = self.ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
                self.model.update({gravity_part_name: self.main_model_part.GetSubModelPart(gravity_part_name)})

        ## Remeshing processes construction
        if (self.ProjectParameters.Has("initial_remeshing_process") == True):
            remeshing_processes = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["initial_remeshing_process"])
            if (ProjectParameters.Has("list_other_processes") == True):
                remeshing_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])

            ## Remeshing processes initialization
            print("STARTING ADAPTATIVE LOOP")
            if (self.ProjectParameters.Has("adaptative_loop") == True):
                adaptative_loop = ProjectParameters["adaptative_loop"].GetInt()
            else:
                adaptative_loop = 1
            for n in range(adaptative_loop):
                print("ADAPTATIVE INTERATION: ", n + 1)
                for process in reversed(remeshing_processes):
                    process.ExecuteInitialize()
                    process.ExecuteBeforeSolutionLoop()
                    process.ExecuteInitializeSolutionStep()

                    if (self.output_post == True):
                        output_settings = ProjectParameters["output_configuration"]
                        gid_output_initial = GiDOutputProcess(self.solver.GetComputingModelPart(),
                                                            self.problem_name+"_"+str(n+1),
                                                            output_settings)
                        gid_output_initial.ExecuteInitialize()
                        gid_output_initial.ExecuteBeforeSolutionLoop()
                        gid_output_initial.ExecuteInitializeSolutionStep()
                        gid_output_initial.ExecuteFinalizeSolutionStep()
                        gid_output_initial.PrintOutput()
                        gid_output_initial.ExecuteFinalize()

        # Obtain the list of the processes to be applied
        if self.problem_type  == "fluid":
            self.list_of_processes = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["initial_conditions_process_list"] )
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )
        elif self.problem_type  == "solid":
            self.list_of_processes = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["constraints_process_list"])
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["loads_process_list"])
        if (self.ProjectParameters.Has("list_other_processes") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])
        if (self.ProjectParameters.Has("json_check_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["json_check_process"])
        if (self.ProjectParameters.Has("json_output_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["json_output_process"])
        if (self.ProjectParameters.Has("compare_two_files_check_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["compare_two_files_check_process"])
        if (self.ProjectParameters.Has("recursive_remeshing_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(self.ProjectParameters["recursive_remeshing_process"])

        for process in self.list_of_processes:
            process.ExecuteInitialize()

        # ### START SOLUTION ####

        self.computing_model_part = self.solver.GetComputingModelPart()
        self.root_model_part = self.computing_model_part.GetRootModelPart()

        if (self.output_post == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

    def Solve(self):
        if self.solve_problem == True:
            for process in self.list_of_processes:
                process.ExecuteBeforeSolutionLoop()

            # #Stepping and time settings (get from process info or solving info)
            # Delta time
            delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
            # Start time
            time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
            # End time
            end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
            step = 0

            if self.problem_type  == "fluid":
                init_step = 3
            elif self.problem_type  == "solid":
                init_step = 1

            # Solving the problem (time integration)
            while(time <= end_time):
                time = time + delta_time
                self.solver.AdvanceInTime(time)
                step = step + 1
                self.root_model_part.ProcessInfo[KratosMultiphysics.STEP] = step

                if(step >= init_step):
                    for process in self.list_of_processes:
                        process.ExecuteInitializeSolutionStep()

                    if (self.root_model_part.Is(KratosMultiphysics.MODIFIED)):
                        # WE INITIALIZE THE SOLVER
                        self.solver.Initialize()
                        # WE RECOMPUTE THE PROCESSES AGAIN
                        ## Processes initialization
                        for process in self.list_of_processes:
                            process.ExecuteInitialize()
                        ## Processes before the loop
                        for process in self.list_of_processes:
                            process.ExecuteBeforeSolutionLoop()
                        ## Processes of initialize the solution step
                        for process in self.list_of_processes:
                            process.ExecuteInitializeSolutionStep()

                    if (self.output_post == True):
                        self.gid_output.ExecuteInitializeSolutionStep()

                    self.solver.Clear()
                    self.solver.Solve()

                    if (self.output_post == True):
                        self.gid_output.ExecuteFinalizeSolutionStep()

                    for process in self.list_of_processes:
                        process.ExecuteFinalizeSolutionStep()

                    for process in self.list_of_processes:
                        process.ExecuteBeforeOutputStep()

                    if (self.output_post == True):
                        if self.gid_output.IsOutputStep():
                            self.gid_output.PrintOutput()

                    for process in self.list_of_processes:
                        process.ExecuteAfterOutputStep()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalize()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
