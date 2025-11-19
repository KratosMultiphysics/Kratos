import numpy as np

import KratosMultiphysics

# Auxiliary function to check the parallel type at runtime
#TODO: Delete this once we come up with the final factory-based design
def _CheckIsDistributed():
    if KratosMultiphysics.ParallelEnvironment.HasDataCommunicator("World"):
        world_data_comm = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")
        return world_data_comm.IsDistributed()
    else:
        return False
# If required, import parallel applications and modules
if _CheckIsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    import KratosMultiphysics.mpi.distributed_import_model_part_utility as distributed_import_model_part_utility
# Importing factories
if _CheckIsDistributed():
    import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as linear_solver_factory
else:
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    import KratosMultiphysics.base_convergence_criteria_factory as convergence_criteria_factory

import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
from KratosMultiphysics.ConvectionDiffusionApplication import topology_optimization_pde_filter_solver

from KratosMultiphysics import assign_vector_by_direction_process

from KratosMultiphysics import DataCommunicator

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTopologyOptimizationPdeFilterProcess(Model, settings["Parameters"])


class ApplyTopologyOptimizationPdeFilterProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings, radius, main_model_part, optimization_model_part, optimization_domain_nodes_mask, nodes_ids_global_to_local_partition_dictionary):
        KratosMultiphysics.Process.__init__(self)
        self.data_communicator = DataCommunicator.GetDefault()
        self.model = Model
        self.main_model_part = main_model_part
        self.domain_size = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        self._SetNumNodesInElements()
        self.optimization_model_part = optimization_model_part
        self.optimization_domain_nodes_mask = optimization_domain_nodes_mask
        self.nodes_ids_global_to_local_partition_dictionary = nodes_ids_global_to_local_partition_dictionary
        self.settings = KratosMultiphysics.Parameters("""
                                    {
                                    "solver_settings": {
                                            "domain_size": """ + str(self.domain_size) + """ 
                                        }
                                    }
                                    """)
        self.filter_radius = radius
        self.filter_reaction = 1.0
        self.pde_solver = topology_optimization_pde_filter_solver.CreateSolver(self.model, self.settings, self.optimization_model_part, self.num_nodes_elements)
        self.pde_solver_set_up = False

    def _SetNumNodesInElements(self):
        ## Elements
        self.num_nodes_elements = 0
        for el in self.main_model_part.Elements:
            self.num_nodes_elements = len(el.GetNodes())
            break
         
    def Execute(self):
        self._InitializePdeFilterExecution()
        self._SolvePdeFilter()
        self._FinalizePdeFilterExecution()       

    def _SolvePdeFilter(self):
        self._GetSolver().InitializeSolutionStep()
        self._GetSolver().Predict()
        self._SynchronizePdeFilterParametersVariables()
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
        if ((not KratosMultiphysics.KratosGlobals.HasVariable("CONVECTION_DIFFUSION_SETTINGS")) or (self._GetSolver().GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS) is None)):
            init_convection_diffusion_settings = KratosMultiphysics.ConvectionDiffusionSettings()
            self._GetSolver().GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, init_convection_diffusion_settings)
        convection_diffusion_settings = self._GetSolver().GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        convection_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable("PDE_FILTER_RESULT"))
        convection_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("PDE_FILTER_FORCING"))
        convection_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("PDE_FILTER_FLUX"))

    def _FinalizePdeFilterExecution(self):
        self.filtered_value_in_opt_nodes = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self.main_model_part.GetCommunicator().LocalMesh().Nodes, KratosCD.PDE_FILTER_RESULT, 0))[self.optimization_domain_nodes_mask]

    def _SetFilterDiffusionAndReactionCoefficients(self):
        self._GetSolver()._UpdateFilterDiffusionVariable(self.filter_radius)
        self._GetSolver()._UpdateFilterReactionVariable(self.filter_reaction)

    def _UpdatePdeFilterSourceTerm(self, source):
        # source is intended to be defined in the main model part nodes
        for node in self.optimization_model_part.GetCommunicator().LocalMesh().Nodes:
            node.SetSolutionStepValue(KratosCD.PDE_FILTER_FORCING, source[self.nodes_ids_global_to_local_partition_dictionary[node.Id]])

    def _GetLastFilteredValue(self):
        return self.filtered_value_in_opt_nodes
    
    def _SynchronizePdeFilterParametersVariables(self):
        self.main_model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosCD.PDE_FILTER_DIFFUSION)
        self.main_model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosCD.PDE_FILTER_REACTION)
        self.main_model_part.GetCommunicator().SynchronizeVariable(KratosCD.PDE_FILTER_FORCING)



