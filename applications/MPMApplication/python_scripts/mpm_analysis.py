
# Importing Kratos Core, Applications and Dependencies
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.MPMApplication.python_solvers_wrapper_mpm import CreateSolver

class MPMAnalysis(AnalysisStage):
    """
    This class is the main-script of the MPMApplication put in a class
    """

    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initializing the parameters
        solver_settings = project_parameters["solver_settings"]

        # Import parallel modules if needed
        # has to be done before the base-class constuctor is called (in which the solver is constructed)
        # TODO: currently MPI parallelization is not present in MPMApplication
        # TODO: the following import lines will be kept here for future reference
        if (project_parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            warn_msg  = 'Currently MPI parallelization is not present in MPMApplication!'
            KratosMultiphysics.Logger.PrintWarning("MPMAnalysis", warn_msg)
            # import KratosMultiphysics.MetisApplication as MetisApplication
            # import KratosMultiphysics.TrilinosApplication as TrilinosApplication

        super(MPMAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes"""
        list_of_processes = super(MPMAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "gravity"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                info_msg  = "Using the old way to create the processes, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintWarning("MPMAnalysis", info_msg)
                from KratosMultiphysics.process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("grid_output_configuration") or self.project_parameters.Has("body_output_configuration"):
                info_msg  = "Using the old way to create the gid-output for grid, this will be removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                KratosMultiphysics.Logger.PrintInfo("MPMAnalysis", info_msg)
                if self.project_parameters.Has("grid_output_configuration"):
                    grid_gid_output= self._SetUpGiDOutput("grid_output")
                    list_of_processes += [grid_gid_output,]
                if self.project_parameters.Has("body_output_configuration"):
                    mp_gid_output= self._SetUpGiDOutput("body_output")
                    list_of_processes += [mp_gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _SetUpGiDOutput(self, parameter_name):
        '''Initialize a GiD output instance'''
        if self.parallel_type == "OpenMP":
            if parameter_name == "grid_output":
                from KratosMultiphysics.gid_output_process import GiDOutputProcess as OutputProcess
                grid_output_file_name = self.project_parameters["problem_data"]["problem_name"].GetString() + "_Grid"
                gid_output = OutputProcess(self._GetSolver().GetGridModelPart(), grid_output_file_name,
                                    self.project_parameters["grid_output_configuration"])
            elif parameter_name == "body_output":
                from KratosMultiphysics.MPMApplication.mpm_gid_output_process import MPMGiDOutputProcess as OutputProcess
                mp_output_file_name = self.project_parameters["problem_data"]["problem_name"].GetString() + "_Body"
                gid_output = OutputProcess(self._GetSolver().GetComputingModelPart(), mp_output_file_name,
                                    self.project_parameters["body_output_configuration"])
        return gid_output

    def _GetSimulationName(self):
        return "::[MPM Analysis]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 mpm_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 mpm_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = MPMAnalysis(model, parameters)
    simulation.Run()
