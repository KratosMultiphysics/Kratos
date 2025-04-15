import numpy as np

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
from KratosMultiphysics.ConvectionDiffusionApplication import topology_optimization_pde_filter_solver

from KratosMultiphysics import assign_vector_by_direction_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTopologyOptimizationPdeFilterProcess(Model, settings["Parameters"])


class ApplyTopologyOptimizationPdeFilterProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings, radius, main_model_part, optimization_model_part, optimization_domain_nodes_mask):
        KratosMultiphysics.Process.__init__(self)
        self.model = Model
        self.main_model_part = main_model_part
        self.domain_size = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        self.optimization_model_part = optimization_model_part
        self.optimization_domain_nodes_mask = optimization_domain_nodes_mask
        self.settings = KratosMultiphysics.Parameters("""
                                    {
                                    "solver_settings": {
                                            "domain_size": """ + str(self.domain_size) + """ 
                                        }
                                    }
                                    """)
        self.filter_radius = radius
        self.filter_reaction = 1.0
        self.pde_solver = topology_optimization_pde_filter_solver.CreateSolver(self.model, self.settings, self.optimization_model_part)
        self.pde_solver_set_up = False

    def Execute(self):
        self._InitializePdeFilterExecution()
        self._SolvePdeFilter()
        self._FinalizePdeFilterExecution()       

    def _SolvePdeFilter(self):
        self._GetSolver().InitializeSolutionStep()
        self._GetSolver().Predict()
        is_converged = self._GetSolver().SolveSolutionStep()
        self._GetSolver().FinalizeSolutionStep()
    
    def _GetSolver(self):
        return self.pde_solver
    
    def _InitializePdeFilterExecution(self):
        if ( not self.pde_solver_set_up):
            self._SetUpPdeFilterSolver()
            self.pde_solver_set_up = True
        self._PreparePdeFilterConvectionDiffusionSettings()
    
    def _SetUpPdeFilterSolver(self):
        self._GetSolver().AddVariables()
        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()
        self._GetSolver().Initialize()
        self._SetFilterDiffusionAndReactionCoefficients()
    
    def _PreparePdeFilterConvectionDiffusionSettings(self):
        convention_diffusion_settings = self._GetSolver().GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable("PDE_FILTER_RESULT"))
        convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("PDE_FILTER_FORCING"))
        convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("PDE_FILTER_FLUX"))

    def _FinalizePdeFilterExecution(self):
        self.filtered_value_in_opt_nodes = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self.main_model_part.Nodes, KratosCD.PDE_FILTER_RESULT, 0))[self.optimization_domain_nodes_mask]

    def _SetFilterDiffusionAndReactionCoefficients(self):
        self._GetSolver()._UpdateFilterDiffusionVariable(self.filter_radius)
        self._GetSolver()._UpdateFilterReactionVariable(self.filter_reaction)

    def _UpdatePdeFilterSourceTerm(self, source):
        # source is intende to be defined in the main model part nodes
        for node in self.optimization_model_part.Nodes:
            node.SetSolutionStepValue(KratosCD.PDE_FILTER_FORCING, source[node.Id-1])

    def _GetLastFilteredValue(self):
        return self.filtered_value_in_opt_nodes



