# Importing the Kratos Library
import KratosMultiphysics

import numpy as np #import the numpy library
import scipy as sp #import the scipy library

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


# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Importing the base class
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_stationary_solver import ConvectionDiffusionStationarySolver
import KratosMultiphysics.FluidDynamicsApplication

def CreateSolver(model, custom_settings, optimization_model_part):
    solver_settings = custom_settings["solver_settings"]
    return TopologyOptimizationPdeFilterSolver(model, solver_settings, optimization_model_part)

class TopologyOptimizationPdeFilterSolver(ConvectionDiffusionStationarySolver):

    def __init__(self, model, custom_settings, optimization_model_part):
        self._DisableSettingsComputeReactions(custom_settings)
        self._SetSettingsModelPartName(custom_settings)
        super().__init__(model,custom_settings)
        self.base_optimization_model_part = optimization_model_part
        self.solver_imports_model_part = False
        self._DefineElementsAndConditions()
    
    def _DisableSettingsComputeReactions(self, custom_settings):
        if not custom_settings.Has("compute_reactions"):
            custom_settings.AddEmptyValue("compute_reactions").SetBool(False)
        else:
            custom_settings["compute_reactions"].SetBool(False)

    def _SetSettingsModelPartName(self, custom_settings):
        if not custom_settings.Has("model_part_name"):
            custom_settings.AddEmptyValue("model_part_name").SetString("PdeFilterModelPart")
        else:
            custom_settings["model_part_name"].SetString("PdeFilterModelPart")

    def _DefineElementsAndConditions(self):
        self.element_name = "TopologyOptimizationPdeFilterElement"
        self.condition_name = "ThermalFace"
        self.element_integrates_in_time = True

    def _get_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")

        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)

        ## Elements
        for el in self.main_model_part.Elements:
            num_nodes_elements = len(el.GetNodes())
            break
        name_string = f"{self.element_name}{domain_size}D{num_nodes_elements}N"
        self.settings["element_replace_settings"]["element_name"].SetString(name_string)
        ## Conditions
        num_nodes_conditions = domain_size
        name_string = f"{self.condition_name}{domain_size}D{num_nodes_conditions}N"
        self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]

    def AddVariables(self):
        #  PDE FILTER VARIABLES
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_RESULT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_FORCING)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_DIFFUSION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_REACTION)
        # self._DefineSettings()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Topology Optimization PDE Filter solver variables added correctly.")
        
    def _DefineSettings(self):
        super().AddVariables()
    
    def _SetTimeSchemeBufferSize(self):
        self.settings["time_scheme"].SetString("bdf2")
        self.min_buffer_size = 1

    def ImportModelPart(self):
        # Here the optimization model part is cloned to be pde filter model part so that the nodes are shared
        element_name, condition_name = self.__GetElementAndConditionNames()
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(
            self.base_optimization_model_part,
            self.main_model_part,
            element_name,
            condition_name)

    def PrepareModelPart(self):
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCD.PDE_FILTER_DIFFUSION, 0.0, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCD.PDE_FILTER_REACTION , 0.0, self.main_model_part.Nodes)
        
        if not self.is_restarted():
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self._get_element_condition_replace_settings()).Execute()
            self._set_and_fill_buffer()

        if (self.settings["echo_level"].GetInt() > 0):
            KratosMultiphysics.Logger.PrintInfo(self.model)

        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]::", "ModelPart prepared for Solver.")

    def _UpdateFilterDiffusionVariable(self, pdf_filter_radius):
        mp = self.GetComputingModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosCD.PDE_FILTER_DIFFUSION, pdf_filter_radius**2)
        
    def _UpdateFilterReactionVariable(self, pdf_filter_reaction):
        mp = self.GetComputingModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosCD.PDE_FILTER_REACTION, pdf_filter_reaction)
    
    def AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append("PDE_FILTER_RESULT")
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Topology Optimization PDE Filter solver DOFs added correctly.")

    def _CreateScheme(self):
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = 1.0
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 0.0
        # As the (no) time integration is managed by the element, we set a "fake" scheme to perform the solution update
        if not _CheckIsDistributed():
            solution_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        else:
            solution_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticScheme()

        return solution_scheme
    
    def SolveSolutionStep(self):
        # Call the base fluid solver to solve current time step
        print("--|PDE_FILTER| ---> Solve Solution Step...")
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        print("--|PDE_FILTER| ---> Step Solved!")
        return is_converged
    
    def AdvanceInTime(self, current_time):
        return super().AdvanceInTime(current_time)
    
    def GetMainModelPart(self):
        return self.main_model_part
    
    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()

    def __GetElementAndConditionNames(self):
        ''' Auxiliary function to get the element and condition names for the connectivity preserve modeler call
        This function returns the element and condition names from the domain size and number of nodes.
        Note that throughout all the substitution process a unique element type and condition is assumed.
        Also note that the connectivity preserve modeler call will create standard base elements as these are to
        be substituted by the corresponding ones in the PrepareModelPart call of the transport solver.
        '''
        ## Get and check domain size
        model_part = self.base_optimization_model_part
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")
        ## Elements
        ## Get the number of nodes from the fluid mesh elements (if there are no elements simplicial are assumed)
        num_nodes_elements = 0
        if (len(model_part.Elements) > 0):
            for elem in model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break
        num_nodes_elements = model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
        if not num_nodes_elements:
            num_nodes_elements = domain_size + 1

        element_name = f"Element{domain_size}D{num_nodes_elements}N"
        ## Conditions
        ## Get the number of nodes from the fluid mesh conditions (if there are no elements simplicial are assumed)
        num_nodes_conditions = 0
        if (len(model_part.Conditions) > 0):
            for cond in model_part.Conditions:
                num_nodes_conditions = len(cond.GetNodes())
                break
        num_nodes_conditions = model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
        if not num_nodes_conditions:
            num_nodes_conditions = domain_size
        aux_condition_name = "LineCondition" if domain_size == 2 else "SurfaceCondition"
        condition_name = f"{aux_condition_name}{domain_size}D{num_nodes_conditions}N"
        return element_name, condition_name

        

