from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSM

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

# Other imports
import sys

# Import the base structural analysis
from structural_mechanics_analysis import StructuralMechanicsAnalysis as BaseClass

class ContactStructuralMechanicsAnalysis(BaseClass):
    """
    This class is the main-script of the ContactStructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):

        # Construct the base analysis.
        super(ContactStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(ContactStructuralMechanicsAnalysis, self).Initialize()
        self.solver.SetEchoLevel(self.echo_level)

    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

        # To avoid many prints
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Solver construction
        import python_solvers_wrapper_contact_structural
        self.solver = python_solvers_wrapper_contact_structural.CreateSolver(self.model, self.project_parameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        self.solver.ImportModelPart() # TODO move to global instance

    def _SetUpListOfProcesses(self):
        """ Set up the list of processes """

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
        if (self.project_parameters.Has("contact_process_list") is True):
            self.list_of_processes += factory.ConstructListOfProcesses(self.project_parameters["contact_process_list"])

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

    def _GetSimulationName(self):
        return "::[KCSM Simulation]:: "

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

    ContactStructuralMechanicsAnalysis(project_parameters_file_name).Run()
