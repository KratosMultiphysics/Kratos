# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication.filtering.python_solvers_wrapper_implicit_filters as implicit_filter_solvers
import KratosMultiphysics.OptimizationApplication as KOA

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class ImplicitFiltering(AnalysisStage):
    """
    This class is the main-script of the implicit filtering analysis
    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases fail with a meaningful error message
        solver_settings = project_parameters["solver_settings"]
        if not solver_settings.Has("time_stepping"):
            raise Exception("Using the old way to pass the time_step, this was removed!")

        if not solver_settings.Has("domain_size"):
            raise Exception("Using the old way to pass the domain_size, this was removed!")

        if not solver_settings.Has("model_part_name"):
            raise Exception("Using the old way to pass the model_part_name, this was removed!")

        if not solver_settings.Has("echo_level"): # this is done to remain backwards-compatible
            raise Exception('"solver_settings" does not have "echo_level", please add it!')

        super().__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        return implicit_filter_solvers.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It was removed, now throwing errors to avoid silent failures
        """
        list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["boundary_conditions_process_list"]
            for process_name in processes_block_names:
                if self.project_parameters.Has(process_name):
                    raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                raise Exception("Using the old way to create the gid-output, this was removed!")
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _GetSimulationName(self):
        return "::[Implicit Helmholtz/Sobolev Filtering]:: "
    
    def FilterField(self, unfiltered_field: KM.ContainerExpression.NodalNonHistoricalExpression) -> KM.ContainerExpression.NodalNonHistoricalExpression:
        unfiltered_field.Read(KOA.HELMHOLTZ_SCALAR_SOURCE)
        unfiltered_field.Evaluate()
        self.Run()
