# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Import applications modules
from KratosMultiphysics.FluidDynamicsApplication import trilinos_fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import trilinos_transport_topology_optimization_solver

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

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

from KratosMultiphysics import DataCommunicator

def CreateSolver(main_model_part, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return TrilinosFluidTransportTopologyOptimizationSolver(main_model_part, solver_settings, is_adjoint_solver=isAdjointSolver)

class TrilinosFluidTransportTopologyOptimizationSolver(PythonSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermallyCoupled",
            "domain_size" : -1,
            "echo_level": 0,
            "fluid_solver_settings": {
                "solver_type": "monolithic",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                }
            },
            "transport_solver_settings": {
                "solver_type": "transient",
                "analysis_type": "linear",
                "model_import_settings"              : {
                    "input_type"     : "mdpa",
                    "input_filename" : "unknown_name"
                },
                "material_import_settings": {
                    "materials_filename": "TransportMaterials.json"
                }
            }
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings, is_adjoint_solver):
        # Initialize MPI utilities
        self._InitializeModelPartImporter()
        self.InitializeDataCommunicator()
        ## Cal base class constructor
        super().__init__(model, custom_settings)
        ## Define if the adjoint problem
        self._DefineAdjointSolver(is_adjoint_solver)
        ## Get domain size
        self._SetDomainSize()
        ## Create subdomain solvers
        self._CreateFluidAndTransportSolvers()

    def InitializeDataCommunicator(self):
        self.data_communicator = DataCommunicator.GetDefault()
        self.nodes_ids_global_to_local_partition_dictionary = {}

    def _InitializeModelPartImporter(self):
        self.distributed_model_part_importer = None

    def _SetDomainSize(self):
        self.domain_size = self.settings["domain_size"].GetInt()
        if not self.settings["fluid_solver_settings"]["solver_settings"].Has("domain_size"):
            self.settings["fluid_solver_settings"]["solver_settings"].AddEmptyValue("domain_size")
        self.settings["fluid_solver_settings"]["solver_settings"]["domain_size"].SetInt(self.domain_size)
        if not self.settings["transport_solver_settings"]["solver_settings"].Has("domain_size"):
            self.settings["transport_solver_settings"]["solver_settings"].AddEmptyValue("domain_size")
        self.settings["transport_solver_settings"]["solver_settings"]["domain_size"].SetInt(self.domain_size)

    def _CreateFluidAndTransportSolvers(self):
        self._CreateFluidSolver()
        self._CreateTransportSolver()

    def _CreateFluidSolver(self):
        self.fluid_solver = trilinos_fluid_topology_optimization_solver.CreateSolver(self.model, self.settings["fluid_solver_settings"], isAdjointSolver=self.IsAdjoint())

    def _CreateTransportSolver(self):
        self.transport_solver = trilinos_transport_topology_optimization_solver.CreateSolver(self.model,self.settings["transport_solver_settings"], isAdjointSolver=self.IsAdjoint()) 
    
    def _DefineAdjointSolver(self, is_adjoint_solver):
        self.is_adjoint = is_adjoint_solver
        
    def AddVariables(self):
        # Import the fluid and transport solver variables. Then merge them to have them in both fluid and transport solvers.
        self.fluid_solver.AddVariables()
        self.transport_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.transport_solver.main_model_part)

    def _DefineTransportProperties(self, decay, convection_coefficient):
        self._GetTransportSolver()._DefineProperties(decay, convection_coefficient)

    def ImportModelPart(self, model_parts=None, physics_solver_distributed_model_part_importer=None):
        if (self.IsPhysics()):
            # FLUID
            self.fluid_solver.ImportModelPart()
            #TRANSPORT
            # Save the convection diffusion settings
            convection_diffusion_settings = self.transport_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
            # Here the fluid model part is cloned to be transport model part so that the nodes are shared
            transport_element_name, transport_condition_name = self.__GetElementAndConditionNames(physics="transport")
            modeler = KratosMultiphysics.ConnectivityPreserveModeler()
            modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.transport_solver.main_model_part, transport_element_name, transport_condition_name)
            # Set the saved convection diffusion settings to the new transport model part
            self.transport_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
            if _CheckIsDistributed():
                self.comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
                ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.transport_solver.main_model_part.GetRootModelPart(), self.comm)
                ParallelFillCommunicator.Execute()
                self.transport_solver.distributed_model_part_importer = self.fluid_solver.distributed_model_part_importer
                self.distributed_model_part_importer = self.fluid_solver.distributed_model_part_importer
        else:
            # FLUID
            fluid_mp = model_parts[0]
            modeler = KratosMultiphysics.ConnectivityPreserveModeler()
            element_name, condition_name = self.__GetElementAndConditionNames(physics="fluid")
            modeler.GenerateModelPart(fluid_mp, self.fluid_solver.main_model_part, element_name, condition_name)
            if _CheckIsDistributed():
                self.comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
                ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.fluid_solver.main_model_part.GetRootModelPart(), self.comm)
                ParallelFillCommunicator.Execute()
                self.fluid_solver.distributed_model_part_importer = physics_solver_distributed_model_part_importer
            # TRANSPORT
            transport_mp = model_parts[1]
            self.fluid_solver.main_model_part = fluid_mp
            self.transport_solver.main_model_part = transport_mp
            if _CheckIsDistributed():
                self.transport_solver.distributed_model_part_importer = physics_solver_distributed_model_part_importer
            # GENERAL
            if _CheckIsDistributed():
                self.distributed_model_part_importer = physics_solver_distributed_model_part_importer

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        default_settings = self.GetDefaultParameters()
        self.settings.ValidateAndAssignDefaults(default_settings)
    
    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        self.transport_solver.PrepareModelPart()

    def _SetAnalysisTimeCoefficient(self):
        self.fluid_solver._SetAnalysisTimeCoefficient()
        self.transport_solver._SetAnalysisTimeCoefficient()

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.transport_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()
    
    def GetMainModelPart(self):
        return self.fluid_solver.GetMainModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_transport = self.transport_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_transport)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.transport_solver.Initialize()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.transport_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.transport_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.transport_solver).SetEchoLevel(level)
        (self.transport_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        if (self.IsPhysics()):
            self.transport_solver.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_T_STEP] += 1
        else:
            self.transport_solver.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_ADJ_T_STEP] += 1
        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.transport_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.transport_solver.Predict()

    def SolveSolutionStep(self):
        if (self.IsPhysics()):
            self.MpiPrint("--|PHY| --|NS| ---> Top. Opt. solution: Solve Fluid Physics", min_echo=0)
            fluid_is_converged = self.fluid_solver.SolveSolutionStep()
            self.MpiPrint("--|PHY| --|T| ---> Top. Opt. solution: Solve Transport Physics", min_echo=0)
            self._SaveVelocityForTransportSolution()
            transport_is_converged = self.transport_solver.SolveSolutionStep()
        elif (self.IsAdjoint()):
            self.MpiPrint("--|PHY| --|ADJ-T| ---> Top. Opt. solution: Solve Transport Adjoint", min_echo=0)
            transport_is_converged = self.transport_solver.SolveSolutionStep()
            self.MpiPrint("--|PHY| --|ADJ_NS| ---> Top. Opt. solution: Solve Fluid Adjoint", min_echo=0)
            self._SaveAdjointTransportScalarForAdjointFluidSolution()
            fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("WRONG TOPOLOGY OPTIMIZATION STAGE", "Running 'SolveSolutionStep' during the wrong topopology optimization stage.")
        return (fluid_is_converged and transport_is_converged)
    
    def _SaveVelocityForTransportSolution(self):
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                KratosMultiphysics.VELOCITY,
                KratosCD.CONVECTION_VELOCITY,
                self._GetFluidSolver().GetMainModelPart(),
                self._GetTransportSolver().GetMainModelPart(),
                0)
        
    def _SaveAdjointTransportScalarForAdjointFluidSolution(self):
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                KratosMultiphysics.TEMPERATURE_ADJ,
                KratosMultiphysics.FUNCTIONAL_DERIVATIVE_TRANSPORT_SCALAR_ADJ,
                self._GetTransportSolver().GetMainModelPart(),
                self._GetFluidSolver().GetMainModelPart(),
                0)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.transport_solver.FinalizeSolutionStep()

    def __GetElementAndConditionNames(self, physics):
        if (physics == "fluid"):
            base_element_name = self.fluid_solver.element_name
            base_condition_name = self.fluid_solver.condition_name
        elif (physics == "transport"):
            base_element_name = self.transport_solver.element_name
            base_condition_name = self.transport_solver.condition_name
        model_part = self.GetMainModelPart()
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
        if not num_nodes_elements:
            num_nodes_elements = domain_size + 1
        element_name = f"{base_element_name}{domain_size}D{num_nodes_elements}N"
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
        condition_name = f"{base_condition_name}{domain_size}D{num_nodes_conditions}N" 
        return element_name, condition_name
    
    def _GetLocalMeshNodes(self):
        return self._GetFluidSolver()._GetLocalMeshNodes()
    
    def _GetFluidSolver(self):
        return self.fluid_solver
    
    def _GetTransportSolver(self):
        return self.transport_solver
    
    def IsAdjoint(self):
        return self.is_adjoint
    
    def IsPhysics(self):
        return not self.IsAdjoint()
    
    def _SetTopologyOptimizationStage(self, top_opt_stage):
        self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TOP_OPT_PROBLEM_STAGE, top_opt_stage)
        self._SetFluidTopologyOptimizationStage(top_opt_stage)
        self._SetTransportTopologyOptimizationStage(top_opt_stage)

    def _GetTopologyOptimizationStage(self):
        return self.GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.TOP_OPT_PROBLEM_STAGE)    
    
    def _SetFluidTopologyOptimizationStage(self, top_opt_stage):
        self._GetFluidSolver()._SetTopologyOptimizationStage(top_opt_stage)

    def _GetFluidTopologyOptimizationStage(self):
        self._GetFluidSolver()._GetTopologyOptimizationStage()

    def _SetTransportTopologyOptimizationStage(self, top_opt_stage):
        self._GetTransportSolver()._SetTopologyOptimizationStage(top_opt_stage)

    def _GetFluidTopologyOptimizationStage(self):
        self._GetTransportSolver()._GetTopologyOptimizationStage()

    def _CheckMaterialProperties(self):
        self._GetFluidSolver()._CheckMaterialProperties()
        self._GetTransportSolver()._CheckMaterialProperties()

    def _UpdateResistanceVariable(self, resistance):
        self._GetFluidSolver()._UpdateResistanceVariable(resistance)

    def _UpdateConductivityVariable(self, conductivity):
        self._GetTransportSolver()._UpdateConductivityVariable(conductivity)

    def _UpdateDecayVariable(self, decay):
        self._GetTransportSolver()._UpdateDecayVariable(decay)

    def _UpdateConvectionCoefficientVariable(self, convection_coefficient):
        self._GetTransportSolver()._UpdateConvectionCoefficientVariable(convection_coefficient)

    def _UpdateTransportSourceVariable(self, transport_source):
        self._GetTransportSolver()._UpdateTransportSourceVariable(transport_source)

    def _UpdateConvectionVelocityVariable(self, convection_velocity):
        self._GetFluidSolver()._UpdateConvectionVelocityVariable(convection_velocity)
        self._GetTransportSolver()._UpdateConvectionVelocityVariable(convection_velocity)

    def _UpdateFunctionalDerivativeVelocityVariable(self, convection_velocity):
        self._GetFluidSolver()._UpdateFunctionalDerivativeVelocityVariable(convection_velocity)
        self._GetTransportSolver()._UpdateFunctionalDerivativeVelocityVariable(convection_velocity)

    def _UpdateFunctionalDerivativeTransportScalarVariable(self, transport_scalar):
        self._GetTransportSolver()._UpdateFunctionalDerivativeTransportScalarVariable(transport_scalar)

    def _UpdateFunctionalDerivativeTransportScalarAdjointVariable(self, transport_scalar_adj):
        for node in self._GetLocalMeshNodes():
            nodal_functional_derivative_transport_scalar_adj = transport_scalar_adj[self.nodes_ids_global_to_local_partition_dictionary[node.Id]]
            node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE_TRANSPORT_SCALAR_ADJ, nodal_functional_derivative_transport_scalar_adj)

    def SetNodesIdsGlobalToLocalDictionary(self, dict):
        self.nodes_ids_global_to_local_partition_dictionary = dict
        self.fluid_solver.SetNodesIdsGlobalToLocalDictionary(dict)
        self.transport_solver.SetNodesIdsGlobalToLocalDictionary(dict)

    def IsNodesIdsGlobalToLocalDictionaryEmpty(self):
        return self.nodes_ids_global_to_local_partition_dictionary == {}
    
    def _IsUnsteady(self):
        return self.fluid_solver._IsUnsteady()
    
    def _IsSteady(self):
        return not self._IsUnsteady()
    
    def _GetTimeSteppingSettings(self):
        return self.fluid_solver._GetTimeSteppingSettings()

    def _SetStartingTime(self, start_time):
        self.fluid_solver.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TIME, start_time)
        self.transport_solver.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TIME, start_time)

    def _SetTimeStep(self, time_step):
        self.fluid_solver._GetTimeSteppingSettings()["time_step"].SetDouble(time_step)
        self.transport_solver._GetTimeSteppingSettings()["time_step"].SetDouble(time_step)

    def MpiCheck(self, text="before solving step", rank=-1):
            if (rank == -1): # print for all ranks
                self.MpiPrint("--|" + str(self.data_communicator.Rank()) + "| Checkpoint reached: " + text)
            elif (self.data_communicator.Rank() == rank): # print only for a specific rank
                self.MpiPrint("--|" + str(rank) + "| Checkpoint reached: " + text)

    def MpiBarrier(self):
        self.data_communicator.Barrier()

    def MpiPrint(self, text_to_print="", rank=0, set_barrier=False, min_echo=1):
        if self.settings["echo_level"].GetInt() >= min_echo:
            if (not _CheckIsDistributed()):
                print(text_to_print)
            else:
                if (set_barrier):
                    self.MpiBarrier()
                if (self.MpiRunOnlyRank(rank)):
                    print(text_to_print)
                if (set_barrier):
                    self.MpiBarrier()    
        else:
            pass 

    def MpiRunOnlyRank(self, rank=0):
        """
        Returns: True if the simulation is not distributed or if it is running on a specified data_communicator rank
        """
        if (not _CheckIsDistributed()):
            return True
        elif (self.data_communicator.Rank() == rank):
            return True
        else:
            return False

    
