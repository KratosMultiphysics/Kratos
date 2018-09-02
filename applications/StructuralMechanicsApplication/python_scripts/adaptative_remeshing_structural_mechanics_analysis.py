from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA

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
try:
    import KratosMultiphysics.MeshingApplication as MA
    KratosMultiphysics.Logger.PrintInfo("MeshingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshingApplication", "not imported")

# Other imports
import sys

# Import the base structural analysis
from structural_mechanics_analysis import StructuralMechanicsAnalysis as BaseClass

class AdaptativeRemeshingStructuralMechanicsAnalysis(BaseClass):
    """
    This class is the main-script of the StructuralMechanicsApplication when using adaptative remeshing put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):

        # Construct the base analysis.
        default_params = KM.Parameters("""
        {
            "max_iteration" : 1,
            "analysis_type" : "linear"
        }
        """)
        if (project_parameters["solver_settings"].Has("max_iteration") is True):
            self.non_linear_iterations = project_parameters["solver_settings"]["max_iteration"].GetInt()
        else:
            self.non_linear_iterations = 10
            project_parameters["solver_settings"].AddValue("max_iteration", default_params["max_iteration"])
        if (project_parameters["solver_settings"].Has("analysis_type") is True):
            project_parameters["solver_settings"]["analysis_type"].SetString("linear")
        else:
            project_parameters["solver_settings"].AddValue("analysis_type", default_params["analysis_type"])
        if (project_parameters.Has("recursive_remeshing_process") is True):
            self.process_remesh = True
        else:
            self.process_remesh = False
        super(AdaptativeRemeshingStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(AdaptativeRemeshingStructuralMechanicsAnalysis, self).Initialize()
        if (self.process_remesh is False):
            convergence_criteria = self._GetSolver().get_convergence_criterion()
            convergence_criteria.Initialize(self._GetSolver().GetComputingModelPart())
        # Ensuring to have conditions on the BC before remesh
        computing_model_part = self._GetSolver().GetComputingModelPart()
        # We need to detect the conditions in the boundary conditions
        if (self.project_parameters.Has("constraints_process_list") is True):
            constraints_process_list = self.project_parameters["constraints_process_list"]
            list_model_parts = []
            for i in range(0,constraints_process_list.size()):
                item = constraints_process_list[i]
                list_model_parts.append(item["Parameters"]["model_part_name"].GetString())
        skin_detection_parameters = KM.Parameters("""
        {
            "list_model_parts_to_assign_conditions" : []
        }
        """)
        for name_mp in list_model_parts:
            skin_detection_parameters["list_model_parts_to_assign_conditions"].Append(name_mp)
        if (computing_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2):
            detect_skin = KM.SkinDetectionProcess2D(computing_model_part, skin_detection_parameters)
        else:
            detect_skin = KM.SkinDetectionProcess3D(computing_model_part, skin_detection_parameters)
        detect_skin.Execute()
        self._GetSolver().SetEchoLevel(self.echo_level)

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        # If we remesh using a process
        if (self.process_remesh is True):
            while self.time < self.end_time:
                self.time = self._GetSolver().AdvanceInTime(self.time)
                if (self.main_model_part.Is(KratosMultiphysics.MODIFIED) is True):
                    # WE INITIALIZE THE SOLVER
                    self._GetSolver().Initialize()
                    # WE RECOMPUTE THE PROCESSES AGAIN
                    ## Processes initialization
                    for process in self._list_of_processes:
                        process.ExecuteInitialize()
                    ## Processes before the loop
                    for process in self._list_of_processes:
                        process.ExecuteBeforeSolutionLoop()
                self.InitializeSolutionStep()
                self._GetSolver().Predict()
                self._GetSolver().SolveSolutionStep()
                self.FinalizeSolutionStep()
                self.OutputSolutionStep()
        else: # Remeshing adaptively
            computing_model_part = self._GetSolver().GetComputingModelPart()
            metric_process = self._GetSolver().get_metric_process()
            remeshing_process = self._GetSolver().get_remeshing_process()
            convergence_criteria = self._GetSolver().get_convergence_criterion()
            builder_and_solver = self._GetSolver().get_builder_and_solver()
            mechanical_solution_strategy = self._GetSolver().get_mechanical_solution_strategy()

            while self.time < self.end_time:
                self.time = self._GetSolver().AdvanceInTime(self.time)
                non_linear_iteration = 1
                while non_linear_iteration <= self.non_linear_iterations:
                    if (computing_model_part.Is(KratosMultiphysics.MODIFIED) is True):
                        self._GetSolver().Clear()
                        # WE RECOMPUTE THE PROCESSES AGAIN
                        # Processes initialization
                        for process in self._list_of_processes:
                            process.ExecuteInitialize()
                        ## Processes before the loop
                        for process in self._list_of_processes:
                            process.ExecuteBeforeSolutionLoop()
                        ## Processes of initialize the solution step
                        for process in self._list_of_processes:
                            process.ExecuteInitializeSolutionStep()
                    if (non_linear_iteration == 1 or computing_model_part.Is(KratosMultiphysics.MODIFIED) is True):
                        self.InitializeSolutionStep()
                        self._GetSolver().Predict()
                        computing_model_part.Set(KratosMultiphysics.MODIFIED, False)
                        self.is_printing_rank = False
                    computing_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, non_linear_iteration)
                    is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self._GetSolver().SolveSolutionStep()
                    is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self.FinalizeSolutionStep()
                    if (is_converged):
                        self.is_printing_rank = True
                        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Adaptative strategy converged in ", non_linear_iteration, "iterations" )
                        break
                    elif (non_linear_iteration == self.non_linear_iterations):
                        self.is_printing_rank = True
                        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Adaptative strategy not converged after ", non_linear_iteration, "iterations" )
                        break
                    else:
                        metric_process.Execute()
                        remeshing_process.Execute()
                        computing_model_part.Set(KratosMultiphysics.MODIFIED, True)
                        non_linear_iteration += 1
                self.OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

        # To avoid many prints
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Solver construction
        import python_solvers_wrapper_adaptative_remeshing_structural
        return python_solvers_wrapper_adaptative_remeshing_structural.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(AdaptativeRemeshingStructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["recursive_remeshing_process"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KratosMultiphysics.Logger.PrintInfo("AdaptativeRemeshingStructuralMechanicsAnalysis", "Using the old way to create the processes, this will be removed!")
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            pass # Already added
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 contact_structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 contact_structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    AdaptativeRemeshingStructuralMechanicsAnalysis(project_parameters_file_name).Run()
