# Kratos imports
import KratosMultiphysics
from KratosMultiphysics.FreeSurfaceApplication.edgebased_levelset_solver import EdgeBasedLevelSetSolver
from KratosMultiphysics.analysis_stage import AnalysisStage


class FreeSurfaceAnalysis(AnalysisStage):
    '''Main script for free surface application.'''

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        # Base class constructor defines:
        # - self.model
        # - self.project_parameters
        # - self.echo_level
        # - self.parallel_type
        super().__init__(model,parameters)

    def Initialize(self):
        super().Initialize()

        model_part_name = self.project_parameters["problem_data"]["problem_name"].GetString()

        # Initialize remaining nodes
        density = self.project_parameters["solver_settings"]["density"].GetDouble()
        body_force = self.project_parameters["solver_settings"]["body_force"].GetVector()
        small_value = 1e-4
        active_node_count = 0
        for node in self.model.GetModelPart(model_part_name).Nodes:
            # Initialize DISTANCE
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                active_node_count += 1
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -small_value)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, small_value)

            # Make sure no node has null porosity and diameter
            if node.GetSolutionStepValue(KratosMultiphysics.POROSITY) == 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.POROSITY, 1.0)
            if node.GetSolutionStepValue(KratosMultiphysics.DIAMETER) == 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DIAMETER, 1.0)

            porosity = node.GetSolutionStepValue(KratosMultiphysics.POROSITY)
            node.SetSolutionStepValue(KratosMultiphysics.PRESS_PROJ, density * body_force * porosity)

        if not active_node_count:
            raise RuntimeError("At least 1 node must be initialized with a negative DISTANCE")

        self._GetSolver()._Redistance()

    def _CreateSolver(self) -> EdgeBasedLevelSetSolver:
        return EdgeBasedLevelSetSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "Free Surface Analysis"