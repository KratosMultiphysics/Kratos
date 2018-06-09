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

class AdaptativeStructuralMechanicsAnalysis(BaseClass):
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
        super(AdaptativeStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(AdaptativeStructuralMechanicsAnalysis, self).Initialize()
        if (self.process_remesh is False):
            convergence_criteria = self.solver.get_convergence_criterion()
            convergence_criteria.Initialize(self.solver.GetComputingModelPart())
        # Ensuring to have conditions on the BC before remesh
        computing_model_part = self.solver.GetComputingModelPart()
        if (computing_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2):
            detect_skin = KratosMultiphysics.SkinDetectionProcess2D(computing_model_part)
        else:
            detect_skin = KratosMultiphysics.SkinDetectionProcess3D(computing_model_part)
        detect_skin.Execute()
        self.solver.SetEchoLevel(self.echo_level)

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        # If we remesh using a process
        if (self.process_remesh is True):
            while self.time < self.end_time:
                self.time = self.solver.AdvanceInTime(self.time)
                if (self.main_model_part.Is(KratosMultiphysics.MODIFIED) is True):
                    # WE INITIALIZE THE SOLVER
                    self.solver.Initialize()
                    # WE RECOMPUTE THE PROCESSES AGAIN
                    ## Processes initialization
                    for process in self.list_of_processes:
                        process.ExecuteInitialize()
                    ## Processes before the loop
                    for process in self.list_of_processes:
                        process.ExecuteBeforeSolutionLoop()
                self.InitializeSolutionStep()
                self.solver.Predict()
                self.solver.SolveSolutionStep()
                self.FinalizeSolutionStep()
                self.OutputSolutionStep()
        else: # Remeshing adaptively
            computing_model_part = self.solver.GetComputingModelPart()
            remeshing_process = self.solver.get_remeshing_process()
            convergence_criteria = self.solver.get_convergence_criterion()
            builder_and_solver = self.solver.get_builder_and_solver()
            mechanical_solution_strategy = self.solver.get_mechanical_solution_strategy()

            while self.time < self.end_time:
                self.time = self.solver.AdvanceInTime(self.time)
                non_linear_iteration = 1
                while non_linear_iteration <= self.non_linear_iterations:
                    if (computing_model_part.Is(KratosMultiphysics.MODIFIED) is True):
                        # Set again all  GiD  I/O
                        self._SetUpGiDOutput()
                        if self.have_output:
                            self.list_of_processes[-1] = self.output
                        # WE RECOMPUTE THE PROCESSES AGAIN
                        # Processes initialization
                        for process in self.list_of_processes:
                            process.ExecuteInitialize()
                        ## Processes before the loop
                        for process in self.list_of_processes:
                            process.ExecuteBeforeSolutionLoop()
                        ## Processes of initialize the solution step
                        for process in self.list_of_processes:
                            process.ExecuteInitializeSolutionStep()
                    if (non_linear_iteration == 1 or computing_model_part.Is(KratosMultiphysics.MODIFIED) is True):
                        self.InitializeSolutionStep()
                        self.solver.Predict()
                        computing_model_part.Set(KratosMultiphysics.MODIFIED, False)
                        self.is_printing_rank = False
                    computing_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, non_linear_iteration)
                    is_converged = convergence_criteria.PreCriteria(computing_model_part, builder_and_solver.GetDofSet(), mechanical_solution_strategy.GetSystemMatrix(), mechanical_solution_strategy.GetSolutionVector(), mechanical_solution_strategy.GetSystemVector())
                    self.solver.SolveSolutionStep()
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
                        remeshing_process.Execute()
                        computing_model_part.Set(KratosMultiphysics.MODIFIED, True)
                        non_linear_iteration += 1
                self.OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

        # To avoid many prints
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Solver construction
        import python_solvers_wrapper_adaptative_structural
        self.solver = python_solvers_wrapper_adaptative_structural.CreateSolver(self.model, self.project_parameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        self.solver.ImportModelPart() # TODO move to global instance

    def _SetUpListOfProcesses(self):
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        self.list_of_processes = factory.ConstructListOfProcesses(self.project_parameters["constraints_process_list"])
        self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["loads_process_list"])
        if (self.project_parameters.Has("list_other_processes") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["list_other_processes"])
        if (self.project_parameters.Has("json_output_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_output_process"])
        # Processes for tests
        if (self.project_parameters.Has("json_check_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["json_check_process"])
        if (self.project_parameters.Has("check_analytic_results_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["check_analytic_results_process"])
        if (self.project_parameters.Has("recursive_remeshing_process") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["recursive_remeshing_process"])

        #TODO this should be generic
        # initialize GiD  I/O
        self._SetUpGiDOutput()
        if self.have_output:
            self.list_of_processes += [self.output,]

        if self.is_printing_rank and self.echo_level > 1:
            count = 0
            for process in self.list_of_processes:
                count += 1
                # KratosMultiphysics.Logger.PrintInfo("Process " + str(count), process) # FIXME
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

    AdaptativeStructuralMechanicsAnalysis(project_parameters_file_name).Run()
