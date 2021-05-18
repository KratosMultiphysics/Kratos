# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if not CheckIfApplicationsAvailable("CompressiblePotentialFlowApplication"):
    raise ImportError("The CompressiblePotentialFlowApplication is not available!")

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_forces_on_nodes_process import ComputeForcesOnNodesProcess
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_lift_process import ComputeLiftProcess

def Create(settings, model, solver_name):
    return PotentialFlowWrapper(settings, model, solver_name)

class PotentialFlowWrapper(kratos_base_wrapper.KratosBaseWrapper):
    def _CreateAnalysisStage(self):
        return PotentialFlowAnalysis(self.model, self.project_parameters)

    def Predict(self):
        pass

    def Initialize(self):

        super().Initialize()

        sub_project_parameters = self.project_parameters["processes"]["boundary_conditions_process_list"]

        for i in range(sub_project_parameters.size()):
            if sub_project_parameters[i]["python_module"].GetString() == "define_wake_process_2d":
                self.wake_process = DefineWakeProcess2D(self.model, sub_project_parameters[i]["Parameters"])
                if not hasattr(self, "wake_process"):
                    raise Exception("potential flow requires specification of a process for the wake (currently specifically using 'define_wake_process_2d')")

            if sub_project_parameters[i]["python_module"].GetString() == "compute_forces_on_nodes_process":
                self.conversion_process = ComputeForcesOnNodesProcess(self.model, sub_project_parameters[i]["Parameters"])
            if sub_project_parameters[i]["python_module"].GetString() == "compute_lift_process":
                self.lift_process = ComputeLiftProcess(self.model, sub_project_parameters[i]["Parameters"])

    def SolveSolutionStep(self):
        self.wake_process.ExecuteInitialize()

        ## the next two lines are needed in order to add Wake DoFs to the new Wake Elements Nodes
        ## and delete the ones that are no longer in the Wake Region.
        self._analysis_stage._GetSolver().Clear()
        self._analysis_stage._GetSolver().InitializeSolutionStep()

        super().SolveSolutionStep()

        self.lift_process.ExecuteFinalizeSolutionStep()
        self.conversion_process.ExecuteFinalizeSolutionStep()
