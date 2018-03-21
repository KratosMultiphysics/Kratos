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
import structural_mechanics_analysis

class ContactStructuralMechanicsAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):
    """
    This class is the main-script of the ContactStructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, project_parameters, external_model_part=None):
        if (type(project_parameters) == str): # a file name is provided
            with open(project_parameters,'r') as parameter_file:
                self.ProjectParameters = KM.Parameters(parameter_file.read())
        elif (type(project_parameters) == KM.Parameters): # a Parameters object is provided
            self.ProjectParameters = project_parameters
        else:
            raise Exception("Input is expected to be provided as a Kratos Parameters object or a file name")
        self.__CreateSolver(external_model_part)

    #### Internal functions ####
    def __CreateSolver(self, external_model_part=None):
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
        import python_solvers_wrapper_structural
        self.solver = python_solvers_wrapper_structural.CreateSolver(self.main_model_part, self.ProjectParameters)

        ## Adds the necessary variables to the model_part only if they don't exist
        self.solver.AddVariables()

        if not self.using_external_model_part:
            ## Read the model - note that SetBufferSize is done here
            self.solver.ReadModelPart() # TODO move to global instance

    def __ExecuteInitialize(self):
        """ Initializing the Analysis """

        ## ModelPart is being prepared to be used by the solver
        self.solver.PrepareModelPartForSolver()

        ## Adds the Dofs if they don't exist
        self.solver.AddDofs()

        ## Creation of the Kratos model (build sub_model_parts or submeshes)
        self.structure_model = KM.Model()
        self.structure_model.AddModelPart(self.main_model_part)

        ## Print model_part and properties
        if self.is_printing_rank and self.echo_level > 1:
            KM.Logger.PrintInfo("ModelPart", self.main_model_part)
            for properties in self.main_model_part.Properties:
                KM.Logger.PrintInfo("Property " + str(properties.Id), properties)

        ## Processes construction
        import process_factory
        self.list_of_processes = process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["constraints_process_list"])
        self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["loads_process_list"])
        if (self.ProjectParameters.Has("list_other_processes") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])
        if (self.ProjectParameters.Has("contact_process_list") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["contact_process_list"])
        if (self.ProjectParameters.Has("json_output_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["json_output_process"])
        # Processes for tests
        if (self.ProjectParameters.Has("json_check_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["json_check_process"])
        if (self.ProjectParameters.Has("check_analytic_results_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["check_analytic_results_process"])

        if self.is_printing_rank and self.echo_level > 1:
            count = 0
            for process in self.list_of_processes:
                count += 1
                # KM.Logger.PrintInfo("Process " + str(count), process) # FIXME

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        ## Add the processes to the solver
        self.solver.AddProcessesList(self.list_of_processes)
        if (self.output_post == True):
            self.solver.AddPostProcess(self.gid_output)

        ## Solver initialization
        self.solver.Initialize()

    def __GetSimulationName(self):
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
