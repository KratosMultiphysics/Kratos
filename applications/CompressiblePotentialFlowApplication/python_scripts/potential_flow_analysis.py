# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class PotentialFlowAnalysis(AnalysisStage):
    '''Main script for potential flow simulations.'''

    def _CreateSolver(self):
        if self.project_parameters["solver_settings"]["solver_type"].GetString()=="potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver
            return potential_flow_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])
        elif self.project_parameters["solver_settings"]["solver_type"].GetString()=="ale_potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.ale_potential_flow_solver as ale_potential_flow_solver
            return ale_potential_flow_solver.CreateSolver(self.model, self.project_parameters["solver_settings"], self.parallel_type)
        elif self.project_parameters["solver_settings"]["solver_type"].GetString()=="adjoint_potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_adjoint_solver as adjoint_solver
            return adjoint_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])
        else:
            raise Exception("Solver type '"+str(self.project_parameters["solver_settings"]["solver_type"].GetString())+"' not added. Please specify an available solver")
