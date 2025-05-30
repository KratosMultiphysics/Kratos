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

from KratosMultiphysics import DataCommunicator

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Importing the base class
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_transient_solver import ConvectionDiffusionTransientSolver
import KratosMultiphysics.FluidDynamicsApplication

def CreateSolver(model, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return TransportTopologyOptimizationSolverMpi(model, solver_settings, is_adjoint_solver=isAdjointSolver)

class TransportTopologyOptimizationSolverMpi(ConvectionDiffusionTransientSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "linear_solver_settings": {
                "solver_type": "amgcl",
                "smoother_type":"ilu0",
                "krylov_type":"gmres",
                "coarsening_type":"aggregation", 
                "max_iteration": 5000,
                "tolerance": 1e-9,
                "scaling": false
                }
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings, is_adjoint_solver):
        self._InitializeModelPartImporter()
        self.InitializeDataCommunicator()
        super().__init__(model,custom_settings)
        self._DefineAdjointSolver(is_adjoint_solver)
        self._DefineElementsAndConditions()
        # self._InitializePhysicsParameters()
        
        print_str = "Construction of TransportTopologyOptimizationSolver "
        if self.IsAdjoint():
            print_str += "for Adjoint problem "
        print_str +=  "finished."
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, print_str)

    def InitializeDataCommunicator(self):
        self.data_communicator = DataCommunicator.GetDefault()
        self.nodes_ids_global_to_local_partition_dictionary = {}

    def _DefineAdjointSolver(self, isAdjointSolver):
        self.is_adjoint = isAdjointSolver

    def _DefineElementsAndConditions(self):
        self.element_name = "TransportTopologyOptimizationElement"
        self.condition_name = "ThermalFace"
        self.element_integrates_in_time = True

    def _InitializeModelPartImporter(self):
        self.distributed_model_part_importer = None

    def _get_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")

        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)

        ## Elements
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            for elem in self.main_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break
        num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
        if not num_nodes_elements:
            num_nodes_elements = domain_size + 1
        name_string = f"{self.element_name}{domain_size}D{num_nodes_elements}N"
        self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        ## Conditions
        num_nodes_conditions = 0
        if (len(self.main_model_part.Conditions) > 0):
            for cond in self.main_model_part.Conditions:
                num_nodes_conditions = len(cond.GetNodes())
                break
        num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
        if not num_nodes_conditions:
            num_nodes_conditions = domain_size
        name_string = f"{self.condition_name}{domain_size}D{num_nodes_conditions}N"
        self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]

    def AddVariables(self): 
        # Add Transport Topology Optimization Variables 
        # DESIGN PARAMETER
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DESIGN_PARAMETER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DESIGN_PARAMETER_GRADIENT) 
        #  PHYSICS
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_GRADIENT)
        # ADJOINT
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.HEAT_FLUX_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.FACE_HEAT_FLUX_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_ADJ_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.OPTIMIZATION_TEMPERATURE)
        # GENERAL
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        # EVENTUAL PDE FILTER VARIABLES
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_RESULT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_FORCING)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_DIFFUSION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_REACTION)
        # SETTINGS
        self._DefineSettings()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid Transport Topology Optimization T and T_ADJ solver variables added correctly.")
        
    def _DefineSettings(self):
        super().AddVariables()
    
    def _CheckMaterialProperties(self):
        for node in self._GetLocalMeshNodes():
            self.MpiPrint("--|--> Conductivity: " + str(node.GetValue(KratosMultiphysics.CONDUCTIVITY)))
            self.MpiPrint("--|--> Decay: " + str(node.GetValue(KratosCD.DECAY)))
            self.MpiPrint("--|--> Convection Coefficient: " + str(node.GetValue(KratosMultiphysics.CONVECTION_COEFFICIENT)))
            break
    
    def _SetTimeSchemeBufferSize(self):
        self.settings["time_scheme"].SetString("bdf2")
        self.min_buffer_size = 1

    def PrepareModelPart(self):
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.CONDUCTIVITY, 0.0, self._GetLocalMeshNodes())
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCD.DECAY, 0.0, self._GetLocalMeshNodes())
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.CONVECTION_COEFFICIENT, 0.0, self._GetLocalMeshNodes())
        if not self.is_restarted():
            # Import material properties
            materials_imported = self.import_materials()
            if materials_imported:
                KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were successfully imported.")
            else:
                KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were not imported.")
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()
            self._set_and_fill_buffer()

        # Create the MPI communicators
        if _CheckIsDistributed():
            self.distributed_model_part_importer.CreateCommunicators()

        if (self.settings["echo_level"].GetInt() > 0):
            KratosMultiphysics.Logger.PrintInfo(self.model)
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]::", "ModelPart prepared for Solver.")

    def _UpdateConductivityVariable(self, conductivity):
        self.is_conductivity_updated = True
        if isinstance(conductivity, (int, float)): 
            for node in self._GetLocalMeshNodes():
             node.SetValue(KratosMultiphysics.CONDUCTIVITY, conductivity)
        elif isinstance(conductivity, (np.ndarray, list)): 
            for node in self._GetLocalMeshNodes():
                node.SetValue(KratosMultiphysics.CONDUCTIVITY, conductivity[self.nodes_ids_global_to_local_partition_dictionary[node.Id]])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateConductivityVariable' : {type(conductivity)}")
        
    def _UpdateDecayVariable(self, decay):
        self.is_decay_updated = True
        if isinstance(decay, (int, float)): 
            for node in self._GetLocalMeshNodes():
             node.SetValue(KratosCD.DECAY, decay)
        elif isinstance(decay, (np.ndarray, list)): 
            for node in self._GetLocalMeshNodes():
                node.SetValue(KratosCD.DECAY, decay[self.nodes_ids_global_to_local_partition_dictionary[node.Id]])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateDecayVariable' : {type(decay)}")

    def _UpdateConvectionCoefficientVariable(self, convection_coefficient):
        self.is_convection_coefficient_updated = True
        if isinstance(convection_coefficient, (int, float)): 
            for node in self._GetLocalMeshNodes():
             node.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, convection_coefficient)
        elif isinstance(convection_coefficient, (np.ndarray, list)): 
            for node in self._GetLocalMeshNodes():
                node.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, convection_coefficient[self.nodes_ids_global_to_local_partition_dictionary[node.Id]])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateResistanceVariable' : {type(convection_coefficient)}")
        
    def _UpdateTransportSourceVariable(self, transport_source):
        self.is_transport_source_updated = True
        if isinstance(transport_source, (int, float)): 
            for node in self._GetLocalMeshNodes():
             node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, transport_source)
        elif isinstance(transport_source, (np.ndarray, list)): 
            for node in self._GetLocalMeshNodes():
                node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, transport_source[self.nodes_ids_global_to_local_partition_dictionary[node.Id]])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateDecayVariable' : {type(transport_source)}")

    def _SetConvectiveVelocity(self, is_constant_velocity, vel):
        self.is_constant_velocity = is_constant_velocity
        self.constant_velocity = vel
        if (self.is_constant_velocity):
            self._SetConstantConvectiveVelocity(self.constant_velocity)
        else:
            self._SetNonConstantConvectiveVelocity()
    
    def _SetConstantConvectiveVelocity(self, velocity):
        for node in self._GetLocalMeshNodes():
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, KratosMultiphysics.Vector(velocity))

    def _SetNonConstantConvectiveVelocity(self):
        pass
    
    def AddDofs(self):
        dofs_and_reactions_to_add = []
        if (self.IsPhysics()): # T
            dofs_and_reactions_to_add.append("TEMPERATURE")
            dofs_str = "T"
        else: # ADJ T
            dofs_and_reactions_to_add.append("TEMPERATURE_ADJ")
            dofs_str = "ADJ T"
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Transport Topology Optimization " + dofs_str + " solver DOFs added correctly.")

    def _CreateScheme(self):
        # Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = self.settings["transient_parameters"]["theta"].GetDouble()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.settings["transient_parameters"]["dynamic_tau"].GetDouble()

        # As the time integration is managed by the element, we set a "fake" scheme to perform the solution update
        if not _CheckIsDistributed():
            transport_top_opt_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        else:
            transport_top_opt_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticScheme()
        return transport_top_opt_scheme
    
    def SolveSolutionStep(self):
        problem_physiscs = self._GetTopologyOptimizationStage()
        if (problem_physiscs == 1):
            problem_phase_str = "TRANS"
        elif (problem_physiscs == 2):
            problem_phase_str = "ADJ-T"
        else: 
            problem_phase_str = "ERROR|"
        if self.IsNodesIdsGlobalToLocalDictionaryEmpty():
            raise RuntimeError("Executing 'SolveSolutionStep' of FluidTopologyOptimizationSolverMpi' with self.nodes_ids_global_to_local_partition_dictionary == { }")
        self.PrintPhysicsParametersUpdateStatus(problem_phase_str)
        self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " solution step...")
        # Call the base fluid solver to solve current time step
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        self.MpiPrint("--|" + problem_phase_str + "| ---> Step Solved!")
        return is_converged
    
    def AdvanceInTime(self, current_time):
        if (not self.IsAdjoint()): # T
            self.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_T_STEP] += 1
            new_time =  super().AdvanceInTime(current_time)
        else: #ADJ
            self.MpiPrint("[WARNING] The Adjoint problem should go backward buth this has not yet been implemented. Since for now it is steady, we advance in time even if it is unnecessary.")
            dt = self.ComputeDeltaTime()
            # new_time = current_time - dt
            new_time = current_time + dt
            self.main_model_part.CloneTimeStep(new_time)
            # self.MpiPrint("\nASK HOW TO HANDLE THIS!!!\n")
            # self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            self.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_ADJ_T_STEP] += 1
        return new_time
    
    def IsAdjoint(self):
        return self.is_adjoint
    
    def IsPhysics(self):
        return not self.IsAdjoint()
    
    def GetMainModelPart(self):
        return self.main_model_part
    
    def _SetTopologyOptimizationStage(self, top_opt_stage):
        self.GetComputingModelPart().ProcessInfo.SetValue(KratosCD.TRANSPORT_TOP_OPT_PROBLEM_STAGE, top_opt_stage)

    def _GetTopologyOptimizationStage(self):
        return self.GetComputingModelPart().ProcessInfo.GetValue(KratosCD.TRANSPORT_TOP_OPT_PROBLEM_STAGE)

    def ImportModelPart(self, model_parts=None, physics_solver_distributed_model_part_importer=None):
        if (self.IsPhysics()):
            if not _CheckIsDistributed():
                self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])
            else:
                self.distributed_model_part_importer = distributed_import_model_part_utility.DistributedImportModelPartUtility(
                    self.main_model_part,
                    self.settings)
                self.distributed_model_part_importer.ImportModelPart()
        else:
            transport_mp = model_parts[0]
            self.main_model_part = transport_mp
            if _CheckIsDistributed():
                self.distributed_model_part_importer = physics_solver_distributed_model_part_importer
    
    def InitializeSolutionStep(self):
        self.is_conductivity_updated = False
        self.is_decay_updated = False
        self.is_convection_coefficient_updated = False
        self.is_transport_source_updated = True
        self._GetSolutionStrategy().InitializeSolutionStep()

    def PrintPhysicsParametersUpdateStatus(self, problem_phase_str):
        if (not self.is_conductivity_updated):
            self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating CONDUCTIVITY variable")
        if (not self.is_decay_updated):
            self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating DECAY variable")
        if (not self.is_convection_coefficient_updated):
            self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating CONVECTION_COEFFICIENT variable")
        if (not self.is_transport_source_updated):
            self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating TRANSPORT SOURCE variable")

    def MpiCheck(self, text="before solving step", rank=-1):
            if (rank == -1): # print for all ranks
                self.MpiPrint("--|" + str(self.data_communicator.Rank()) + "| Checkpoint reached: " + text)
            elif (self.data_communicator.Rank() == rank): # print only for a specific rank
                self.MpiPrint("--|" + str(rank) + "| Checkpoint reached: " + text)

    def MpiBarrier(self):
        self.data_communicator.Barrier()

    def MpiPrint(self, text_to_print="", rank=0, set_barrier=False):
        if (not _CheckIsDistributed()):
            print(text_to_print)
        else:
            if (set_barrier):
                self.MpiBarrier()
            if (self.MpiRunOnlyRank(rank)):
                print(text_to_print)
            if (set_barrier):
                self.MpiBarrier()  

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

    def SetNodesIdsGlobalToLocalDictionary(self, dict):
        self.nodes_ids_global_to_local_partition_dictionary = dict

    def IsNodesIdsGlobalToLocalDictionaryEmpty(self):
        return self.nodes_ids_global_to_local_partition_dictionary == {}
    
    def _GetLocalMeshNodes(self, mp = None):
        if mp is None:
            return self.GetMainModelPart().GetCommunicator().LocalMesh().Nodes
        else:
            return mp.GetCommunicator().LocalMesh().Nodes

                
            