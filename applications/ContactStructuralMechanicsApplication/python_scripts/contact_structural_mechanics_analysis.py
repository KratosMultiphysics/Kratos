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
    def __init__(self, project_parameters, external_model_part=None):

        # Construct the base analysis.
        super().__init__(project_parameters, external_model_part)

    #### Internal functions ####
    def _CreateSolver(self, external_model_part=None):
        """ Create the Solver (and create and import the ModelPart if it is not passed from outside) """
        if external_model_part != None:
            # This is a temporary solution until the importing of the ModelPart
            # is removed from the solver (needed e.g. for Optimization)
            if (type(external_model_part) != KM.ModelPart):
                raise Exception("Input is expected to be provided as a Kratos ModelPart object")
            self.using_external_model_part = True
        else:
            self.using_external_model_part = False

        ## Get echo level and parallel type
        self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.ProjectParameters["problem_data"]["parallel_type"].GetString()

        # To avoid many prints # TODO leave this?
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KM.mpi as KratosMPI
            import KM.MetisApplication as MetisApplication
            import KM.TrilinosApplication as TrilinosApplication
            self.is_printing_rank = (KratosMPI.mpi.rank == 0)
        else:
            self.is_printing_rank = True

        ## Structure model part definition
        if self.using_external_model_part:
            self.main_model_part = external_model_part
        else:
            main_model_part_name = self.ProjectParameters["problem_data"]["model_part_name"].GetString()
            self.main_model_part = KM.ModelPart(main_model_part_name)
            self.main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE,
                                                      self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        ## Solver construction
        import python_solvers_wrapper_contact_structural
        self.solver = python_solvers_wrapper_contact_structural.CreateSolver(self.main_model_part, self.ProjectParameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        if not self.using_external_model_part:
            ## Read the model - note that SetBufferSize is done here
            self.solver.ReadModelPart() # TODO move to global instance

    def _ExecuteInitialize(self):
        """ Initializing the Analysis """

        super()._ExecuteInitialize()

        ## Processes construction
        import process_factory
        if (self.ProjectParameters.Has("contact_process_list") is True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["contact_process_list"])

            #if self.is_printing_rank and self.echo_level > 1: # FIXME
                #KM.Logger.PrintInfo("Process " + str(len(self.list_of_processes)), self.list_of_processes[-1])

            self.list_of_processes[-1].ExecuteInitialize()

        ## Add the processes to the solver
        self.solver.AddProcessesList(self.list_of_processes)
        if (self.output_post is True):
            self.solver.AddPostProcess(self.gid_output)

        # Initialize the solver (again)
        self.solver.Initialize()

        # Setting the echo level
        echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()
        self.solver.SetEchoLevel(echo_level)

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
