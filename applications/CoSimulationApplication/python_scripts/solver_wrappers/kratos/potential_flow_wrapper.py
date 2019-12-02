from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if not CheckIfApplicationsAvailable("CompressiblePotentialFlowApplication"):
    raise ImportError("The CompressiblePotentialFlowApplication is not available!")
import KratosMultiphysics.CompressiblePotentialFlowApplication

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_forces_on_nodes_process import ComputeForcesOnNodesProcess
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_lift_process import ComputeLiftProcess
from KratosMultiphysics.gid_output_process import Factory as GiDFactory

def Create(settings, solver_name):
    return PotentialFlowWrapper(settings, solver_name)

class PotentialFlowWrapper(kratos_base_wrapper.KratosBaseWrapper):
    def _CreateAnalysisStage(self):
        return PotentialFlowAnalysis(self.model, self.project_parameters)

    def AdvanceInTime(self, current_time):
        self.time = 0.0
        return self.time
    def Predict(self):
        return

    def Initialize(self):

        super(PotentialFlowWrapper, self).Initialize()

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

        gid_parameters = self.project_parameters["output_processes"]["gid_output"][0]
        gid_parameters["Parameters"]["output_name"].SetString('Fsi_Potential_Flow_Output_with_TimeSteps')
        self.gid_output_process = GiDFactory(gid_parameters,self.model)
        self.gid_output_process.ExecuteInitialize()
        self.gid_output_process.ExecuteBeforeSolutionLoop()
        self.time = 1

    def SolveSolutionStep(self):
        self.wake_process.ExecuteInitialize()
        self.gid_output_process.ExecuteInitializeSolutionStep()
        self._analysis_stage._GetSolver().fluid_solver.solver.Clear()
        self._analysis_stage._GetSolver().fluid_solver.solver.InitializeSolutionStep()
        super(PotentialFlowWrapper, self).SolveSolutionStep()
        self.model["FluidModelPart"].CloneTimeStep(self.time)
        self.model["FluidModelPart"].ProcessInfo[KratosMultiphysics.STEP] += 1

        self.lift_process.ExecuteFinalizeSolutionStep()
        self.conversion_process.ExecuteFinalizeSolutionStep()
        self.gid_output_process.ExecuteFinalizeSolutionStep()
        self.gid_output_process.PrintOutput()
        self.time +=1

    def FinalizeSolutionStep(self):
        super(PotentialFlowWrapper, self).FinalizeSolutionStep()
        self.gid_output_process.ExecuteFinalize()
