from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Other imports
import sys

# Import the base structural analysis
from KratosMultiphysics.ContactStructuralMechanicsApplication.contact_structural_mechanics_analysis import ContactStructuralMechanicsAnalysis as BaseClass

# Import solver wrapper
from KratosMultiphysics.ContactStructuralMechanicsApplication import python_solvers_wrapper_adaptative_remeshing_contact_structural

# Import auxiliar methods
from KratosMultiphysics.ContactStructuralMechanicsApplication.auxiliar_methods_contact_adaptative_remeshing import AuxiliarMethodsContactAdaptiveRemeshing

class AdaptativeRemeshingContactStructuralMechanicsAnalysis(BaseClass):
    """
    This class is the main-script of the ContactStructuralMechanicsApplication when using adaptative remeshing put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):

        # Construct the base analysis.
        default_params = KM.Parameters("""
        {
            "max_iteration"    : 1,
            "analysis_type"    : "linear",
            "contact_settings" :
            {
                "fancy_convergence_criterion" : false
            }
        }
        """)
        if project_parameters["problem_data"].Has("generate_auxiliar_boundary"):
            self.generate_auxiliar_boundary = project_parameters["problem_data"]["generate_auxiliar_boundary"].GetBool()
        else:
            self.generate_auxiliar_boundary = False
        if project_parameters["solver_settings"].Has("max_iteration"):
            self.non_linear_iterations = project_parameters["solver_settings"]["max_iteration"].GetInt()
        else:
            self.non_linear_iterations = 10
            project_parameters["solver_settings"].AddValue("max_iteration", default_params["max_iteration"])
        if not project_parameters["solver_settings"].Has("analysis_type"):
            project_parameters["solver_settings"].AddValue("analysis_type", default_params["analysis_type"])
        if project_parameters["solver_settings"].Has("contact_settings"):
            if project_parameters["solver_settings"]["contact_settings"].Has("fancy_convergence_criterion"):
                project_parameters["solver_settings"]["contact_settings"]["fancy_convergence_criterion"].SetBool(False)
            else:
                project_parameters["solver_settings"]["contact_settings"].AddValue("fancy_convergence_criterion", default_params["contact_settings"]["fancy_convergence_criterion"])
        else:
            project_parameters["solver_settings"].AddValue("contact_settings", default_params["contact_settings"])
        self.process_remesh = ""
        if project_parameters.Has("mesh_adaptivity_processes"):
            self.process_remesh = "remesh_loop"
        if project_parameters.Has("processes"):
            if project_parameters["processes"].Has("mesh_adaptivity_processes"):
                self.process_remesh = "remesh_loop"
        if project_parameters["solver_settings"].Has("convergence_criterion"):
            if project_parameters["solver_settings"]["convergence_criterion"].GetString() == "adaptative_remesh_criteria":
                self.process_remesh = "adaptively"
        if self.process_remesh == "adaptively":
            project_parameters["solver_settings"]["analysis_type"].SetString("linear")
        super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

        # Create utilities
        self.adaptive_utilities = AuxiliarMethodsContactAdaptiveRemeshing(self)

    def Initialize(self):
        """ Initializing the Analysis """
        super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).Initialize()

        # If the boundary is detected automatically
        if self.generate_auxiliar_boundary :
            self.adaptive_utilities.AdaptativeRemeshingDetectBoundary()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        # If we remesh using a process
        if self.process_remesh == "remesh_loop":
            self.adaptive_utilities.AdaptativeRemeshingRunSolutionLoop()
        elif self.process_remesh == "adaptively": # Remeshing adaptively
            self.adaptive_utilities.SPRAdaptativeRemeshingRunSolutionLoop()
        else:
            super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self).RunSolutionLoop()

    def ClearDatabase(self):
        """ This method clears the database in case it is necessary

            Keyword arguments:
            self It signifies an instance of a class.
        """
        solver = self._GetSolver()
        computing_model_part = solver.GetComputingModelPart()
        computing_model_part.ProcessInfo.SetValue(CSMA.ACTIVE_SET_COMPUTED, False)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

        # To avoid many prints
        if self.echo_level == 0:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Solver construction
        return python_solvers_wrapper_adaptative_remeshing_contact_structural.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(AdaptativeRemeshingContactStructuralMechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["mesh_adaptivity_processes"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KM.Logger.PrintWarning("AdaptativeRemeshingContactStructuralMechanicsAnalysis", "Using the old way to create the processes, this will be removed!")
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
        err_msg += '    "python3 adaptative_remeshing_contact_structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 adaptative_remeshing_contact_structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    AdaptativeRemeshingContactStructuralMechanicsAnalysis(project_parameters_file_name).Run()
