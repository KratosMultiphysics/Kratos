# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication.python_solvers_wrapper_mesh_motion as mesh_mothion_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class MeshMovingAnalysis(AnalysisStage):
    """
    This class is the main-script of the MeshMovingApplication put in a class
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
        return mesh_mothion_solvers.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It was removed, now throwing errors to avoid silent failures
        """
        list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["boundary_conditions_process_list", "list_other_processes", "json_output_process",
                "json_check_process", "check_analytic_results_process"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        raise Exception("Using the old way to create the processes, this was removed!")
            else: # Processes are given in the new format
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
        return "::[Mesh Moving Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 mesh_moving_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 mesh_moving_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = MeshMovingAnalysis(model, parameters)
    simulation.Run()
