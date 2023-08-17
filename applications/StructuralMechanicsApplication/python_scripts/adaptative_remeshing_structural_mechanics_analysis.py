# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_adaptative_remeshing_structural

# Import the base structural analysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis as BaseClass

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
        if project_parameters["solver_settings"].Has("max_iteration"):
            self.non_linear_iterations = project_parameters["solver_settings"]["max_iteration"].GetInt()
        else:
            self.non_linear_iterations = 10
            project_parameters["solver_settings"].AddValue("max_iteration", default_params["max_iteration"])
        if project_parameters["solver_settings"].Has("analysis_type"):
            project_parameters["solver_settings"]["analysis_type"].SetString("linear")
        else:
            project_parameters["solver_settings"].AddValue("analysis_type", default_params["analysis_type"])
        self.process_remesh = False
        if project_parameters.Has("mesh_adaptivity_processes"):
            self.process_remesh = True
        if project_parameters.Has("processes"):
            if project_parameters["processes"].Has("mesh_adaptivity_processes"):
                self.process_remesh = True
        super().__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super().Initialize()
        computing_model_part = self._GetSolver().GetComputingModelPart()
        if not self.process_remesh:
            convergence_criteria = self._GetSolver()._GetConvergenceCriterion()
            convergence_criteria.Initialize(computing_model_part)

        # Ensuring to have conditions on the BC before remesh
        is_surface = False
        for elem in computing_model_part.Elements:
            geom = elem.GetGeometry()
            if geom.WorkingSpaceDimension() != geom.LocalSpaceDimension():
                is_surface = True
            break

        if not is_surface:
            # We need to detect the conditions in the boundary conditions
            list_model_parts = []
            if self.project_parameters.Has("constraints_process_list"):
                constraints_process_list = self.project_parameters["constraints_process_list"]
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

            if computing_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
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
        computing_model_part = self._GetSolver().GetComputingModelPart()
        root_model_part = computing_model_part.GetRootModelPart()
        if self.process_remesh:
            while self.time < self.end_time:
                self.time = self._GetSolver().AdvanceInTime(self.time)
                # We reinitialize if remeshed previously
                if root_model_part.Is(KM.MODIFIED):
                    self._ReInitializeSolver()
                self.InitializeSolutionStep()
                # We reinitialize if remeshed on the InitializeSolutionStep
                if root_model_part.Is(KM.MODIFIED):
                    self._ReInitializeSolver()
                    self.InitializeSolutionStep()
                self._GetSolver().Predict()
                self._GetSolver().SolveSolutionStep()
                self.FinalizeSolutionStep()
                self.OutputSolutionStep()
        else: # Remeshing adaptively
            metric_process = self._GetSolver().get_metric_process()
            remeshing_process = self._GetSolver().get_remeshing_process()
            convergence_criteria = self._GetSolver()._GetConvergenceCriterion()
            builder_and_solver = self._GetSolver()._GetBuilderAndSolver()
            mechanical_solution_strategy = self._GetSolver()._GetSolutionStrategy()

            while self.time < self.end_time:
                self.time = self._GetSolver().AdvanceInTime(self.time)
                non_linear_iteration = 1
                while non_linear_iteration <= self.non_linear_iterations:
                    if root_model_part.Is(KM.MODIFIED):
                        self._ReInitializeSolver()
                    if non_linear_iteration == 1 or root_model_part.Is(KM.MODIFIED):
                        self.InitializeSolutionStep()
                        self._GetSolver().Predict()
                        computing_model_part.Set(KM.MODIFIED, False)
                    computing_model_part.ProcessInfo.SetValue(KM.NL_ITERATION_NUMBER, non_linear_iteration)
                    reform_dofs = mechanical_solution_strategy.GetReformDofSetAtEachStepFlag()
                    mechanical_solution_strategy.SetReformDofSetAtEachStepFlag(True)
                    is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self._GetSolver().SolveSolutionStep()
                    is_converged = convergence_criteria.PostCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self.FinalizeSolutionStep()
                    mechanical_solution_strategy.SetReformDofSetAtEachStepFlag(reform_dofs)
                    if is_converged:
                        KM.Logger.PrintInfo(self._GetSimulationName(), "Adaptative strategy converged in ", non_linear_iteration, "iterations" )
                        break
                    elif non_linear_iteration == self.non_linear_iterations:
                        KM.Logger.PrintInfo(self._GetSimulationName(), "Adaptative strategy not converged after ", non_linear_iteration, "iterations" )
                        break
                    else:
                        metric_process.Execute()
                        remeshing_process.Execute()
                        computing_model_part.Set(KM.MODIFIED, True)
                        non_linear_iteration += 1
                self.OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

        # To avoid many prints
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Solver construction
        return python_solvers_wrapper_adaptative_remeshing_structural.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["mesh_adaptivity_processes"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KM.Logger.PrintWarning("AdaptativeRemeshingStructuralMechanicsAnalysis", "Using the old way to create the processes, this will be removed!")
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            pass # Already added
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _ReInitializeSolver(self):
        """ This reinitializes after remesh """
        self._GetSolver().Clear()
        # WE INITIALIZE THE SOLVER
        self._GetSolver().Initialize()
        # WE RECOMPUTE THE PROCESSES AGAIN
        ## Processes initialization
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()
        ## Processes before the loop
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()
        ## Processes of initialize the solution step
        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

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
