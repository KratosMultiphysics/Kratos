from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

# Other imports
import sys

# Import the base structural analysis
from fluid_dynamics_analysis import FluidDynamicsAnalysis as BaseClass

class AdaptativeRemeshingFluidDynamicsAnalysis(BaseClass):
    """
    This class is the main-script of the FluidDynamicsApplication when using adaptative remeshing put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):

        # Construct the base analysis.
        super(AdaptativeRemeshingFluidDynamicsAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """
        super(AdaptativeRemeshingFluidDynamicsAnalysis, self).Initialize()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        # If we remesh using a process
        computing_model_part = self._GetSolver().GetComputingModelPart()
        root_model_part = computing_model_part.GetRootModelPart()

        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            # We reinitialize if remeshed previously
            if root_model_part.Is(KratosMultiphysics.MODIFIED):
                self._ReInitializeSolver()
            self.InitializeSolutionStep()
            # We reinitialize if remeshed on the InitializeSolutionStep
            if root_model_part.Is(KratosMultiphysics.MODIFIED):
                self._ReInitializeSolver()
                self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()


    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(AdaptativeRemeshingFluidDynamicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["recursive_remeshing_process"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KratosMultiphysics.Logger.PrintWarning("AdaptativeRemeshingFluidDynamicsAnalysis", "Using the old way to create the processes, this will be removed!")
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
        err_msg += '    "python3 adaptative_remeshing_fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 adaptative_remeshing_fluid_dynamics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    AdaptativeRemeshingFluidDynamicsAnalysis(project_parameters_file_name).Run()
