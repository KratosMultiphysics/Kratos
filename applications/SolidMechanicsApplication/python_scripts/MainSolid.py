# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import system python modules
import time as timer
import sys
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid


sys.stdout.flush()

class Solution(object):

    def __init__(self, file_parameters="ProjectParameters.json", file_name=None):

        # Time control starts
        print(timer.ctime())

        sys.stdout.flush()

        # Measure process time
        self.t0p = timer.clock()

        # Measure wall time
        self.t0w = timer.time()

        # Import input
        self.ProjectParameters = self._import_project_parameters(file_parameters)

        # Set input file name
        self._set_input_file_name(file_name)

        # Set logger severity level
        self._set_severity_level()

        # Defining the number of threads
        num_threads = self._get_parallel_size()
        if self.ProjectParameters.Has("problem_data"):
            if self.ProjectParameters["problem_data"].Has("threads"):
                num_threads = self.ProjectParameters["problem_data"]["threads"].GetInt()
                self._set_parallel_size(num_threads)

        print(" ")
        print(self._class_prefix()+" [OMP USING", num_threads, "THREADS]")

    def Run(self):

        self.Initialize()

        self.Solve()

        self.Finalize()

    def Initialize(self):

        # Start model
        self.model = self._get_model()

        # Start solver
        self.solver = self._get_solver()

        # Start processes
        self.processes = self._get_processes()

        # Get Solution Step Variables
        self._set_solution_step_variables()

        # Initialize model (read and import)
        self.model.ExecuteInitialize()

        sys.stdout.flush()

        # Set time settings
        self._get_time_settings()

        # Initialize Solver
        self.solver.SetEchoLevel(self.echo_level)
        self.solver.ExecuteInitialize()

        # Import materials
        self.main_model_part = self.model.GetMainModelPart()
        if self._is_not_restarted():
            self._import_materials()

        # Initiliaze processes
        self.processes.ExecuteInitialize()

        # Initialize graphical output (GiD)
        output_model_part = self.model.GetOutputModelPart()
        self.output = self._get_graphical_output(output_model_part)
        self.output.ExecuteInitialize()

        # Adaptive solution
        if self.time_process is not None:
            self.time_process.ExecuteInitialize()

    def Solve(self):

        # Initialize solution Loop
        self.InitializeSolutionLoop()

        # Solving the problem (time integration)
        while self.time < self.end_time:

            self.InitializeSolutionStep()
            if self.SolveSolutionStep():
                self.FinalizeSolutionStep()

            if self.echo_level >= 0:
                sys.stdout.flush()

    def InitializeSolutionLoop(self):

        # Processes to be executed before solution loop
        self.processes.ExecuteBeforeSolutionLoop()

        # First execution before solution loop
        self.model.ExecuteBeforeSolutionLoop()

        # Writing a initial state results file or single file (if no restart)
        if self._is_not_restarted():
            self.output.ExecuteBeforeSolutionLoop()

        # First execution before solution loop
        self.solver.ExecuteBeforeSolutionLoop()

        # Print model_part and properties
        if self.echo_level > 0:
            print("")
            print(self.main_model_part)
            for properties in self.main_model_part.Properties:
                print(properties)

        print(" ")
        print(self._class_prefix()+" Analysis -START- ")

        sys.stdout.flush()

    def PredictTimeStep(self):
        # Predict time step from time integration
        if self.time_process is not None:
            self.time_process.Execute()
            self.delta_time = self.process_info[KratosMultiphysics.DELTA_TIME]
            self.time = self.process_info[KratosMultiphysics.TIME]
            self.step = self.process_info[KratosMultiphysics.STEP]
        else:
            self.delta_time = self.process_info[KratosMultiphysics.DELTA_TIME]
            self.time = self.time + self.delta_time
            self.step = self.step + 1
            self.process_info[KratosMultiphysics.STEP] = self.step
            self.main_model_part.CloneTimeStep(self.time)

        if self.echo_level >= 0:
            print("  [STEP:"+str(self.step)+" TIME:"+"{0:1.{1}f}".format(self.time, 5)+"]")

    def InitializeSolutionStep(self):

        clock_time = self._start_time_measuring()

        # Predict time step
        self.PredictTimeStep()

        # Processes to be executed at the begining of the solution step
        self.processes.ExecuteInitializeSolutionStep()

        # Execution at the begining of the solution step
        self.model.ExecuteInitializeSolutionStep()

        # Execution at the begining of the solution step
        self.output.ExecuteInitializeSolutionStep()

        self._stop_time_measuring(clock_time, "Initialize Step", self.report);

    def SolveSolutionStep(self):

        clock_time = self._start_time_measuring()

        # All steps included (1)(2)(3)
        converged = self.solver.Solve()

        # Step by step (1)
        # self.solver.InitializeSolutionStep()

        # Step by step (2)
        # converged = self.solver.SolveSolutionStep()

        # Step by step (3)
        # self.solver.FinalizeSolutionStep()

        iterations = self.process_info[KratosMultiphysics.NL_ITERATION_NUMBER]
        print("  (-ITER:"+str(iterations)+" CPU:%.2f" % round((timer.clock()-clock_time), 2)+"s-)")

        self._stop_time_measuring(clock_time, "Solve Step", self.report)

        self.process_info[KratosSolid.CONVERGENCE_ACHIEVED] = converged

        return converged

    def FinalizeSolutionStep(self):

        clock_time = self._start_time_measuring()

        # Execution at the end of the solution step
        self.output.ExecuteFinalizeSolutionStep()

        # Processes to be executed at the end of the solution step
        self.processes.ExecuteFinalizeSolutionStep()

        # Execution at the end of the solution step
        self.model.ExecuteFinalizeSolutionStep()

        # Processes to be executed before witting the output
        self.processes.ExecuteBeforeOutputStep()

        # Execution before witting the output
        self.model.ExecuteBeforeOutputStep()

        # Write output results GiD: (frequency writing is controlled internally)
        self._print_output()

        # Processes to be executed after witting the output
        self.processes.ExecuteAfterOutputStep()

        # Execution before witting the output
        self.model.ExecuteAfterOutputStep()

        self._stop_time_measuring(clock_time, "Finalize Step", self.report)


    def Finalize(self):

        # Ending the problem (time integration finished)
        self.output.ExecuteFinalize()

        self.processes.ExecuteFinalize()

        print(self._class_prefix()+" Analysis -END- ")
        print(" ")

        # Check solving information for any problem
        # self.solver.InfoCheck() # InfoCheck not implemented yet.

        # Measure process time
        tfp = timer.clock()

        # Measure wall time
        tfw = timer.time()

        print(self._class_prefix()+" [Elapsed Time = %.2f" % (tfw - self.t0w), "seconds] (%.2f" % (tfp - self.t0p), "seconds of cpu/s time)")
        print(timer.ctime())

    #### Main internal methods ####

    def _print_output(self):
        if self.ProjectParameters.Has("output_configuration"):
            if self.output.IsOutputStep():
                self.output.PrintOutput()
                # Write EigenValues
        if self.process_info.Has(KratosSolid.EIGENVALUE_VECTOR):
            current_vals = [ev for ev in self.main_model_part.ProcessInfo[KratosSolid.EIGENVALUE_VECTOR]]
            print(" EIGENVALUES ", current_vals)

    def _import_project_parameters(self, input_file):
        import input_manager
        self.input_manager = input_manager.InputManager(input_file)
        return self.input_manager.GetProjectParameters()

    def _set_input_file_name(self, file_name):
        if file_name is not None:
            if self.ProjectParameters.Has("problem_data") is False:
                void_parameters = KratosMultiphysics.Parameters("{}")
                self.ProjectParameters.AddValue("problem_data", void_parameters)

            if self.ProjectParameters["problem_data"].Has("problem_name"):
                self.ProjectParameters["problem_data"]["problem_name"].SetString(file_name)
            else:
                self.ProjectParameters["problem_data"].AddEmptyValue("problem_name").SetString(file_name)

            if self.ProjectParameters["model_settings"].Has("input_file_settings") is False:
                void_parameters = KratosMultiphysics.Parameters("{}")
                self.ProjectParameters["model_settings"].AddValue("input_file_settings", void_parameters)

            if self.ProjectParameters["model_settings"]["input_file_settings"].Has("name"):
                self.ProjectParameters["model_settings"]["input_file_settings"]["name"].SetString(file_name)
            else:
                self.ProjectParameters["model_settings"]["input_file_settings"].AddEmptyValue("name").SetString(file_name)

    def _is_not_restarted(self):
        if self.process_info.Has(KratosMultiphysics.IS_RESTARTED):
            if self.process_info[KratosMultiphysics.IS_RESTARTED]:
                return False
            else:
                return True
        else:
            return True

    def _set_solution_step_variables(self):
        solver_variables = self.solver.GetVariables()
        self.model.SetVariables(solver_variables)

        processes_variables = self.processes.GetVariables()
        self.model.SetVariables(processes_variables)

    def _get_model(self):
        import model_manager
        return model_manager.ModelManager(self.ProjectParameters["model_settings"])

    def _get_solver(self):
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        return solver_module.CreateSolver(self.ProjectParameters["solver_settings"]["Parameters"], self.model.GetModel())

    def _get_time_settings(self):

        self.process_info = self.model.GetProcessInfo()

        # Get time parameters
        if self._is_not_restarted():
            if self.ProjectParameters.Has("time_settings"):
                time_settings = self.ProjectParameters["time_settings"]
                time_increment = 1.0
                if self.ProjectParameters["time_settings"].Has("time_step"):
                    time_increment = time_settings["time_step"].GetDouble()
                    self.process_info.SetValue(KratosMultiphysics.DELTA_TIME, time_increment)
                    initial_time = 0.0
                if self.ProjectParameters["time_settings"].Has("start_time"):
                    initial_time = time_settings["start_time"].GetDouble()
                    self.process_info.SetValue(KratosMultiphysics.TIME, initial_time)

        # Set time parameters
        self.step = self.process_info[KratosMultiphysics.STEP]
        self.time = self.process_info[KratosMultiphysics.TIME]

        self.delta_time = self.process_info[KratosMultiphysics.DELTA_TIME]

        self.end_time = self.time + self.delta_time
        self.time_process = None
        if self.ProjectParameters.Has("time_settings"):
            if self.ProjectParameters["time_settings"].Has("end_time"):
                self.end_time = self.ProjectParameters["time_settings"]["end_time"].GetDouble()
            if self.ProjectParameters["time_settings"].Has("time_process"):
                process = self.ProjectParameters["time_settings"]["time_process"]
                process["Parameters"].AddValue("end_time",self.end_time)
                process["Parameters"].AddValue("start_time",self.time)
                process["Parameters"].AddValue("time_step",self.delta_time)
                kratos_module = __import__(process["kratos_module"].GetString())
                python_module = __import__(process["python_module"].GetString())
                self.time_process = python_module.Factory(process, self.GetModel())

    def _import_materials(self):
        # Assign material to model_parts (if Materials.json exists)
        import process_factory

        if(os.path.isfile("Materials.json") or self.input_manager.HasMaterialFile()):

            MaterialParameters = self.input_manager.GetMaterialParameters()

            if MaterialParameters.Has("material_models_list"):

                import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterials

                domain_model = self.model.GetModel()

                assign_materials_processes = process_factory.KratosProcessFactory(domain_model).ConstructListOfProcesses(MaterialParameters["material_models_list"])

                for process in assign_materials_processes:
                    process.Execute()

        elif os.path.isfile("materials.py"):  # legacy

            import constitutive_law_python_utility as constitutive_law_utils

            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.main_model_part, self.process_info[KratosMultiphysics.SPACE_DIMENSION])

            constitutive_law.Initialize()

            print("::[-----Material------]:: Reading file: materials.py ")

        else:
            print("No Materials.json or Materials.py found ")


    def _get_processes(self):
        # Obtain the list of the processes to be applied
        import process_handler

        # get processes parameters
        processes_parameters = self._get_processes_parameters()

        domain_model = self.model.GetModel()
        return process_handler.ProcessHandler(domain_model, processes_parameters)

    def _get_processes_parameters(self):

        processes_parameters = KratosMultiphysics.Parameters("{}")
        processes_parameters.AddEmptyValue("echo_level").SetInt(self.echo_level)
        if self.ProjectParameters.Has("constraints_process_list"):
            processes_parameters.AddValue("constraints_process_list", self.ProjectParameters["constraints_process_list"])
        if self.ProjectParameters.Has("loads_process_list"):
            processes_parameters.AddValue("loads_process_list", self.ProjectParameters["loads_process_list"])
        if self.ProjectParameters.Has("problem_process_list"):
            processes_parameters.AddValue("problem_process_list", self.ProjectParameters["problem_process_list"])
        if self.ProjectParameters.Has("output_process_list"):
            processes_parameters.AddValue("output_process_list", self.ProjectParameters["output_process_list"])
        if self.ProjectParameters.Has("check_process_list"):
            processes_parameters.AddValue("check_process_list", self.ProjectParameters["check_process_list"])

        return processes_parameters

    def _get_graphical_output(self, output_model_part):
        # Output settings start
        if self.ProjectParameters.Has("output_configuration"):
            import gid_output_process
            self.output_settings = self.ProjectParameters["output_configuration"]
            problem_name = "results_output"
            if self.ProjectParameters["problem_data"].Has("problem_name"):
                problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()
            else:
                print(" problem name not supplied -> generic name used : results_output ")
                print(self._class_prefix()+" Output Ready [File: "+problem_name+".*.post.* ]")
            return gid_output_process.GiDOutputProcess(output_model_part, problem_name, self.output_settings)
        else:
            print(self._class_prefix()+" No Output")
            return KratosMultiphysics.Process()

    def _set_severity_level(self):
        # Set echo level
        self.echo_level = 0
        if self.ProjectParameters.Has("problem_data"):
            if self.ProjectParameters["problem_data"].Has("echo_level"):
                self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()

        self.severity = KratosMultiphysics.Logger.Severity.WARNING
        if self.echo_level == 1:
            self.severity = KratosMultiphysics.Logger.Severity.INFO
        elif self.echo_level == 2:
            self.severity = KratosMultiphysics.Logger.Severity.DETAIL
        elif self.echo_level == 3:
            self.severity = KratosMultiphysics.Logger.Severity.DEBUG
        elif self.echo_level == 4:
            self.severity = KratosMultiphysics.Logger.Severity.TRACE

        # Print solving time
        self.report = False
        if self.echo_level > 0:
            self.report = True

        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(self.severity)

    def _set_parallel_size(self, num_threads):
        parallel = KratosMultiphysics.OpenMPUtils()
        parallel.SetNumThreads(int(num_threads))

    def _get_parallel_size(self):
        parallel = KratosMultiphysics.OpenMPUtils()
        return parallel.GetNumThreads()

    def _start_time_measuring(self):
        # Measure process time
        time_ip = timer.clock()
        return time_ip

    def _stop_time_measuring(self, time_ip, process, report):
        # Measure process time
        time_fp = timer.clock()
        if report:
            used_time = time_fp - time_ip
            print(self._class_prefix()+" [ %.2f" % round(used_time, 2), "s", process, " ] ")

    @classmethod
    def _class_prefix(self):
        header = "::[---KSM Simulation--]::"
        return header


if __name__ == "__main__":
    Solution().Run()
