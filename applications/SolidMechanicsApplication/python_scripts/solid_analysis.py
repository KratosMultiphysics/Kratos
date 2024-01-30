
# Import system python modules
import time as timer
import sys
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

sys.stdout.flush()

# Importing the analysis_stage base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class Solution(AnalysisStage):

    def __init__(self, Model, project_parameters="ProjectParameters.json", file_name=None):
        # Time control starts
        print(timer.ctime())

        sys.stdout.flush()

        # Measure process time
        self.t0p = timer.clock()

        # Measure wall time
        self.t0w = timer.time()

        # Import input
        self._project_parameters = self._import_project_parameters(project_parameters)

        # Set input file name
        self._set_input_file_name(file_name)

        # Set logger severity level
        self._set_severity_level()

        # Start model manager
        self._model = self._get_model(Model)

        # Defining the number of threads
        num_threads = self._get_parallel_size()
        if self._project_parameters.Has("problem_data"):
            if self._project_parameters["problem_data"].Has("threads"):
                num_threads = self._project_parameters["problem_data"]["threads"].GetInt()
                self._set_parallel_size(num_threads)

        print(" ")
        print(self._class_prefix()+" [OMP USING", num_threads, "THREADS]")

        # Set solver-model-processes manager initialization (this is done in the Initialize)
        # self._initialize_and_import_project()


    def RunSolutionLoop(self):
        # Solving the problem (time integration)
        while self.time < self.end_time:

            self.InitializeSolutionStep()

            if self.SolveSolutionStep():
                self.FinalizeSolutionStep()
            else:
                self.FinalizeSolutionStep()
                #self.FinalizeNonConvergedStep()

            if self.echo_level >= 0:
                sys.stdout.flush()

    def Initialize(self):
        # Set solver-model-processes manager initialization
        self._initialize_and_import_project()

        # Set time settings
        self._get_time_settings()

        # Initialize Solver
        self._solver.SetEchoLevel(self.echo_level)
        self._solver.ExecuteInitialize()

        # Import materials
        self.main_model_part = self._model.GetMainModelPart()
        if self._is_not_restarted():
            self._import_materials()

        # Initiliaze processes
        self._processes.ExecuteInitialize()

        # Initialize graphical output (GiD)
        output_model_part = self._model.GetOutputModelPart()
        self._output = self._get_graphical_output(output_model_part)
        self._output.ExecuteInitialize()

        # Adaptive solution
        if self.time_process is not None:
            self.time_process.ExecuteInitialize()

        # Initialize solution Loop
        self.InitializeSolutionLoop()

    def InitializeSolutionLoop(self):
        # Processes to be executed before solution loop
        self._processes.ExecuteBeforeSolutionLoop()

        # First execution before solution loop
        self._model.ExecuteBeforeSolutionLoop()

        # Writing a initial state results file or single file (if no restart)
        if self._is_not_restarted():
            print(self._class_prefix()+" Print Initial State ")
            self._output.ExecuteBeforeSolutionLoop()
        else:
            # Write output results GiD: (frequency writing is controlled internally)
            self._output.PrintOutput()
            print(self._class_prefix()+" [ Print Restart State ("+str(self.process_info[KratosMultiphysics.PRINTED_STEP])+") ]")

        # First execution before solution loop
        self._solver.ExecuteBeforeSolutionLoop()
        self.process_info[KratosSolid.CONVERGENCE_ACHIEVED] = True

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
        else:
            self.time += self.process_info[KratosMultiphysics.DELTA_TIME]
            self.process_info[KratosMultiphysics.STEP] += 1
            self.main_model_part.CloneTimeStep(self.time)

        self.time_step = self.process_info[KratosMultiphysics.DELTA_TIME]
        self.time = self.process_info[KratosMultiphysics.TIME]
        self.step = self.process_info[KratosMultiphysics.STEP]

        if self.echo_level >= 0:
            print(" [STEP:"+str(self.step)+" TIME:"+"{:.2e}".format(self.time, 5)+"] "+" (-dt:"+"{:.2e}".format(self.time_step)+"/end:"+"{:.1f}".format(self.end_time)+"-)")

    def InitializeSolutionStep(self):
        #clock_time = self._start_time_measuring()

        # Predict time step
        self.PredictTimeStep()

        # Processes to be executed at the begining of the solution step
        self._processes.ExecuteInitializeSolutionStep()

        # Execution at the begining of the solution step
        self._model.ExecuteInitializeSolutionStep()

        # Execution at the begining of the solution step
        self._output.ExecuteInitializeSolutionStep()

        #self._stop_time_measuring(clock_time, "Initialize Step", self.report);

    def SolveSolutionStep(self):
        clock_time = self._start_time_measuring()

        # All steps included (1)(2)(3)
        converged = self._solver.Solve()

        # Step by step (1)
        # self._solver.InitializeSolutionStep()

        # Step by step (2)
        # converged = self._solver.SolveSolutionStep()

        # Step by step (3)
        # self._solver.FinalizeSolutionStep()

        iterations = self.process_info[KratosMultiphysics.NL_ITERATION_NUMBER]
        print("  (-ITER:"+str(iterations)+" CPU:%.2f" % round((timer.clock()-clock_time), 2)+"s-)")

        #self._stop_time_measuring(clock_time, "Solve Step", self.report)

        self.process_info[KratosSolid.CONVERGENCE_ACHIEVED] = converged

        return converged

    def FinalizeSolutionStep(self):
        #clock_time = self._start_time_measuring()

        # Execution at the end of the solution step
        self._output.ExecuteFinalizeSolutionStep()

        # Processes to be executed at the end of the solution step
        self._processes.ExecuteFinalizeSolutionStep()

        # Execution at the end of the solution step
        self._model.ExecuteFinalizeSolutionStep()

        # Execute output at the end of the solution step
        self.OutputSolutionStep()

        self.non_converged_steps -= 1
        #self._stop_time_measuring(clock_time, "Finalize Step", self.report)

    def FinalizeNonConvergedStep(self):
        if self.process_info[KratosSolid.DIVERGENCE_ACHIEVED]: # force reduction if possible
            self.non_converged_steps = 0

        if (self.non_converged_steps > 0) and (self.time_step < self.maximum_time_step): # continue anyway
            self.non_converged_steps -= 1
            self.process_info[KratosSolid.CONVERGENCE_ACHIEVED] = True
            print(self._class_prefix()+" [ not converged step : preserving time_step ]")
            self._solver.FinalizeSolutionStep()
            self.FinalizeSolutionStep()
        else:
            self.non_converged_steps = self.time_change_delay
            print(self._class_prefix()+" [ not converged step : reducing time_step ]")

    def OutputSolutionStep(self):
        # Processes to be executed before witting the output
        self._processes.ExecuteBeforeOutputStep()

        # Execution before witting the output
        self._model.ExecuteBeforeOutputStep()

        # Write output results GiD: (frequency writing is controlled internally)
        self._print_output()

        # Processes to be executed after witting the output
        self._processes.ExecuteAfterOutputStep()

        # Execution before witting the output
        self._model.ExecuteAfterOutputStep()


    def Finalize(self):
        # Ending the problem (time integration finished)
        self._output.ExecuteFinalize()

        self._processes.ExecuteFinalize()

        print(self._class_prefix()+" Analysis -END- ")
        print(" ")

        # Check solving information for any problem
        # self._solver.InfoCheck() # InfoCheck not implemented yet.

        # Measure process time
        tfp = timer.clock()

        # Measure wall time
        tfw = timer.time()

        print(self._class_prefix()+" [Elapsed Time = %.2f" % (tfw - self.t0w), "seconds] (%.2f" % (tfp - self.t0p), "seconds of cpu/s time)")
        print(timer.ctime())


    #### Main internal methods ####

    def _initialize_and_import_project(self):
        # Start solver
        self._solver = self._get_solver()

        # Start processes
        self._processes = self._get_processes()

        # Get Solution Step Variables
        self._set_solution_step_variables()

        # Initialize model (read and import)
        self._model.ExecuteInitialize()

        sys.stdout.flush()

    def _print_output(self):
        if self._project_parameters.Has("output_configuration"):
            if self._output.IsOutputStep():
                self._output.PrintOutput()
                print(self._class_prefix()+" [ Print Output ("+str(self.process_info[KratosMultiphysics.PRINTED_STEP])+") ]")
                # Write EigenValues
        if self.process_info.Has(KratosSolid.EIGENVALUE_VECTOR):
            current_vals = [ev for ev in self.main_model_part.ProcessInfo[KratosSolid.EIGENVALUE_VECTOR]]
            print(" EIGENVALUES ", current_vals)

    def _import_project_parameters(self, input_file):
        import KratosMultiphysics.SolidMechanicsApplication.input_manager as input_manager
        self.input_manager = input_manager.InputManager(input_file)
        return self.input_manager.GetProjectParameters()

    def _set_input_file_name(self, file_name):
        if file_name is not None:
            if self._project_parameters.Has("problem_data") is False:
                void_parameters = KratosMultiphysics.Parameters("{}")
                self._project_parameters.AddValue("problem_data", void_parameters)

            if self._project_parameters["problem_data"].Has("problem_name"):
                self._project_parameters["problem_data"]["problem_name"].SetString(file_name)
            else:
                self._project_parameters["problem_data"].AddEmptyValue("problem_name").SetString(file_name)

            if self._project_parameters["model_settings"].Has("input_file_settings") is False:
                void_parameters = KratosMultiphysics.Parameters("{}")
                self._project_parameters["model_settings"].AddValue("input_file_settings", void_parameters)

            if self._project_parameters["model_settings"]["input_file_settings"].Has("name"):
                self._project_parameters["model_settings"]["input_file_settings"]["name"].SetString(file_name)
            else:
                self._project_parameters["model_settings"]["input_file_settings"].AddEmptyValue("name").SetString(file_name)

    def _is_restarted(self):
        if self.process_info[KratosMultiphysics.IS_RESTARTED]:
            return True
        else:
            return False

    def _is_not_restarted(self):
        return not self._is_restarted()

    def _set_solution_step_variables(self):
        solver_variables = self._solver.GetVariables()
        self._model.SetVariables(solver_variables)

        processes_variables = self._processes.GetVariables()
        self._model.SetVariables(processes_variables)


    def _get_model(self, Model):
        import KratosMultiphysics.SolidMechanicsApplication.model_manager as model_manager
        return model_manager.ModelManager(Model, self._project_parameters["model_settings"])

    def _get_solver(self):
        #if self._project_parameters["solver_settings"].Has("kratos_module"):
        #    kratos_module = __import__(self._project_parameters["solver_settings"]["kratos_module"].GetString())
        #else:
        #    import KratosMultiphysics.SolversApplication
        import KratosMultiphysics
        if self._project_parameters["solver_settings"].Has("kratos_module"):
            kratos_module_name = self._project_parameters["solver_settings"]["kratos_module"].GetString()
        else:
            kratos_module_name = "KratosMultiphysics.SolidMechanicsApplication"

        solver_module_name = kratos_module_name+"."+self._project_parameters["solver_settings"]["solver_type"].GetString()

        #solver_module = __import__(solver_module_name)
        import importlib
        solver_module = importlib.import_module(solver_module_name)

        return solver_module.CreateSolver(self._project_parameters["solver_settings"]["Parameters"], self._model.GetModel())

    def _get_time_settings(self):
        self.process_info = self._model.GetProcessInfo()

        # Get time parameters
        if self._is_not_restarted():
            if self._project_parameters.Has("time_settings"):
                time_settings = self._project_parameters["time_settings"]
                if self._project_parameters["time_settings"].Has("time_step"):
                    self.process_info.SetValue(KratosMultiphysics.DELTA_TIME, time_settings["time_step"].GetDouble())
                if self._project_parameters["time_settings"].Has("start_time"):
                    self.process_info.SetValue(KratosMultiphysics.TIME, time_settings["start_time"].GetDouble())
        else:
            if self._project_parameters.Has("time_settings"):
                time_settings = self._project_parameters["time_settings"]
                if self._project_parameters["time_settings"].Has("time_step"):
                    self.process_info.SetValue(KratosMultiphysics.DELTA_TIME, time_settings["time_step"].GetDouble())

        # Set time parameters
        self.step = self.process_info[KratosMultiphysics.STEP]
        self.time = self.process_info[KratosMultiphysics.TIME]
        self.time_step = self.process_info[KratosMultiphysics.DELTA_TIME]

        self.end_time = self.time + self.time_step
        self.time_process = None
        self.non_converged_steps = 0
        self.time_change_delay = 4
        self.maximum_time_step = 0
        if self._project_parameters.Has("time_settings"):
            if self._project_parameters["time_settings"].Has("end_time"):
                self.end_time = self._project_parameters["time_settings"]["end_time"].GetDouble()
                self.maximum_time_step = self._project_parameters["time_settings"]["time_step"].GetDouble()
            if self._project_parameters["time_settings"].Has("time_process"):
                process = self._project_parameters["time_settings"]["time_process"]
                if process.Has("steps_update_delay"):
                    self.time_change_delay = process["steps_update_delay"].GetDouble()
                process["Parameters"].AddEmptyValue("end_time").SetDouble(self.end_time)
                process["Parameters"].AddEmptyValue("start_time").SetDouble(self.time)
                process["Parameters"].AddEmptyValue("time_step").SetDouble(self.time_step)
                kratos_module = __import__(process["kratos_module"].GetString())
                python_module = __import__(process["python_module"].GetString())
                self.time_process = python_module.Factory(process, self._model.GetModel())

    def _import_materials(self):
        # Assign material to model_parts (if Materials.json exists)
        import KratosMultiphysics.process_factory

        if(os.path.isfile("Materials.json") or self.input_manager.HasMaterialFile()):

            MaterialParameters = self.input_manager.GetMaterialParameters()

            if MaterialParameters.Has("material_models_list"):

                import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterials

                domain_model = self._model.GetModel()

                assign_materials_processes = KratosMultiphysics.process_factory.KratosProcessFactory(domain_model).ConstructListOfProcesses(MaterialParameters["material_models_list"])

                for process in assign_materials_processes:
                    process.Execute()

                self._model.CleanModel()

        elif os.path.isfile("materials.py"):  # legacy

            import KratosMultiphysics.SolidMechanicsApplication.constitutive_law_python_utility as constitutive_law_utils

            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.main_model_part, self.process_info[KratosMultiphysics.SPACE_DIMENSION])

            constitutive_law.Initialize()

            self._model.CleanModel()

            print("::[-----Material------]:: Reading file: materials.py ")

        else:
            print("No Materials.json or Materials.py found ")


    def _get_processes(self):
        # Obtain the list of the processes to be applied
        import KratosMultiphysics.SolidMechanicsApplication.process_handler as process_handler

        # get processes parameters
        processes_parameters = self._get_processes_parameters()

        domain_model = self._model.GetModel()
        return process_handler.ProcessHandler(domain_model, processes_parameters)

    def _get_processes_parameters(self):

        processes_parameters = KratosMultiphysics.Parameters("{}")
        processes_parameters.AddEmptyValue("echo_level").SetInt(self.echo_level)

        process_lists = ["constraints_process_list","loads_process_list","problem_process_list","output_process_list","check_process_list"]
        for list_name in process_lists:
            if self._project_parameters.Has(list_name):
                processes_parameters.AddValue(list_name, self._project_parameters[list_name])

        return processes_parameters

    def _get_graphical_output(self, output_model_part):

        # Time parameters
        if self._is_restarted():
            import KratosMultiphysics
            self.process_info[KratosMultiphysics.PRINTED_STEP] -= 1

        # Output settings start
        if self._project_parameters.Has("output_configuration"):
            import KratosMultiphysics.gid_output_process
            self._output_settings = self._project_parameters["output_configuration"]
            problem_name = "results_output"
            if self._project_parameters["problem_data"].Has("problem_name"):
                problem_name = self._project_parameters["problem_data"]["problem_name"].GetString()
            else:
                print(" problem name not supplied -> generic name used : results_output ")
                print(self._class_prefix()+" Output Ready [File: "+problem_name+".*.post.* ]")
            return KratosMultiphysics.gid_output_process.GiDOutputProcess(output_model_part, problem_name, self._output_settings)
        else:
            print(self._class_prefix()+" No Output")
            import KratosMultiphysics
            return KratosMultiphysics.Process()


    def _set_severity_level(self):
        # Set echo level
        self.echo_level = 0
        if self._project_parameters.Has("problem_data"):
            if self._project_parameters["problem_data"].Has("echo_level"):
                self.echo_level = self._project_parameters["problem_data"]["echo_level"].GetInt()

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
