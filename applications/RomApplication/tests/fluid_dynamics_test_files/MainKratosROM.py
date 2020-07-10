import KratosMultiphysics
from KratosMultiphysics.RomApplication.fluid_dynamics_analysis_rom import FluidDynamicsAnalysisROM

import sys
import time
import numpy as np

class TestFluidDynamicsROM(FluidDynamicsAnalysisROM):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super().__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()
        self.selected_time_step_solution_container = []

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        this_time_step_results = []
        if np.isclose(self.time,4):
            for node in self._solver.GetComputingModelPart().Nodes:
                this_time_step_results.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
                this_time_step_results.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
                this_time_step_results.append(node.GetSolutionStepValue(KratosMultiphysics.PRESSURE))
            self.selected_time_step_solution_container.append(this_time_step_results)

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

    def EvaluateQuantityOfInterest(self):
        SnapshotMatrix = np.zeros((len(self.selected_time_step_solution_container[0]), len(self.selected_time_step_solution_container)))
        for i in range(len(self.selected_time_step_solution_container)):
            Snapshot_i= np.array(self.selected_time_step_solution_container[i])
            SnapshotMatrix[:,i] = Snapshot_i.transpose()
        return SnapshotMatrix


if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = TestFluidDynamicsROM(model,parameters)
    simulation.Run()