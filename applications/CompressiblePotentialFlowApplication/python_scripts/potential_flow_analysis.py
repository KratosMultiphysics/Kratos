# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp

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
        elif self.project_parameters["solver_settings"]["solver_type"].GetString()=="stochastic_adjoint_potential_flow":
            import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_stochastic_adjoint_solver as adjoint_solver
            return adjoint_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])
        else:
            raise Exception("Solver type '"+str(self.project_parameters["solver_settings"]["solver_type"].GetString())+"' not added. Please specify an available solver")

    def ModifyInitialProperties(self):
        # READING refernce chord from LE and TE position
        if self._GetSolver().main_model_part.HasSubModelPart("LeadingEdgeNode"):
            le = [node for node in self._GetSolver().main_model_part.GetSubModelPart("LeadingEdgeNode").Nodes]
            te = [node for node in self._GetSolver().main_model_part.GetSubModelPart("TrailingEdgeNode").Nodes]
            import math
            self._GetSolver().reference_chord = math.sqrt((le[0].X-te[0].X)**2+(le[0].Y-te[0].Y)**2)
            self._GetSolver().main_model_part.ProcessInfo.SetValue(KCPFApp.REFERENCE_CHORD,self._GetSolver().reference_chord)
        else:
            self._GetSolver().reference_chord = self.project_parameters["solver_settings"]["reference_chord"].GetDouble()

