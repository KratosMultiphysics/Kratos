# Importing the Kratos Library
import KratosMultiphysics
from importlib import import_module

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
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.trilinos_navier_stokes_monolithic_solver import TrilinosNavierStokesMonolithicSolver

def CreateSolver(model, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return TrilinosFluidTopologyOptimizationSolver(model, solver_settings, is_adjoint_solver=isAdjointSolver)

class TrilinosFluidTopologyOptimizationSolver(TrilinosNavierStokesMonolithicSolver):
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
                "max_iteration": 1000,
                "tolerance": 1e-6,
                "scaling": false
                }
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings, is_adjoint_solver):
        self._InitializeModelPartImporter()
        self.InitializeDataCommunicator()
        self.CheckModelPartImportSettings(custom_settings)
        self._DefineAdjointSolver(is_adjoint_solver)
        init_model = self._GetInitModel(model)
        super().__init__(init_model,custom_settings)
        self._DefineElementsAndConditions()
        print_str = "Construction of FluidTopologyOptimizationSolver "
        if self.IsAdjoint():
            print_str += "for Adjoint problem "
        print_str +=  "finished."
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, print_str)

    def InitializeDataCommunicator(self):
        self.data_communicator = DataCommunicator.GetDefault()
        self.nodes_ids_global_to_local_partition_dictionary = {}

    def _InitializeModelPartImporter(self):
        self.distributed_model_part_importer = None

    def CheckModelPartImportSettings(self, settings):
        if (settings["model_import_settings"]["input_type"].GetString() != "mdpa"):
            raise RuntimeError("Executing MPI with the wrong model part import settings")
        
    def _GetInitModel(self, model):
        if self.IsAdjoint():
            return KratosMultiphysics.Model()
        else:
            return model

    def _DefineAdjointSolver(self, isAdjointSolver):
        self.is_adjoint = isAdjointSolver

    def _DefineElementsAndConditions(self):
        if (self.element_name != "FluidTopologyOptimizationElement"):
            self.MpiPrint("[WARNING] " + self.__class__.__name__ + " element_name: \' " + self.element_name + " \' is not compatible with FluidTopologyOptimization. Its value has been reset to default value: \' FluidTopologyOptimizationElement \'")
            self.element_name = "FluidTopologyOptimizationElement"
        self.condition_name = "FluidTopologyOptimizationWallCondition"
        self.element_integrates_in_time = True
        
    def _SetFormulation(self):
        self.element_name = "FluidTopologyOptimizationElement"
        self.condition_name = "FluidTopologyOptimizationWallCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list.append(KratosCFD.RESISTANCE)
        self.non_historical_nodal_properties_variables_list.append(KratosMultiphysics.CONVECTION_COEFFICIENT)
        self.non_historical_nodal_properties_variables_list.append(KratosCFD.CONVECTION_VELOCITY)

    def AddVariables(self):  
        #Add parent class variables
        super().AddVariables()   
        # Add DESIGN PARAMETER
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DESIGN_PARAMETER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DESIGN_PARAMETER_GRADIENT) 
        # Add VELOCITY GRADIENT
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Y_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Z_GRADIENT)
        # Add Adjoint-NS Variables
        # self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ACCELERATION_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.MESH_VELOCITY_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.BODY_FORCE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        # eventual pde filter variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_RESULT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_FORCING)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_DIFFUSION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.PDE_FILTER_REACTION)
        # eventual transport coupling
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_ADJ_GRADIENT)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid Topology Optimization ADJ-NS solver variables added correctly.")          

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bossak":
            KratosMultiphysics.Logger.PrintWarning("[WARNING] " + self.__class__.__name__ + " time scheme_type: \' " + scheme_type + " \' is not compatible with the current implementation of FluidTopologyOptimization. Its value has been reset to default value: \' bdf2 \'")
            self.settings["time_scheme"].SetString("bdf2")
            self.min_buffer_size = 3
        elif scheme_type == "bdf2":
            self.min_buffer_size = 3
        elif scheme_type == "steady":
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
    
    def _SetNodalProperties(self):
        set_density = KratosMultiphysics.DENSITY in self.historical_nodal_properties_variables_list
        set_viscosity = KratosMultiphysics.VISCOSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KratosMultiphysics.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list
        set_resistance = KratosCFD.RESISTANCE in self.non_historical_nodal_properties_variables_list
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                if rho <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            # Get DYNAMIC_VISCOSITY from properties and calculate the kinematic one (VISCOSITY)
            if set_viscosity:
                dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
                if dyn_viscosity <= 0.0:
                    raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
                kin_viscosity = dyn_viscosity / rho
            # Get SOUND_VELOCITY
            if set_sound_velocity:
                if el.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                    sound_velocity = el.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
                else:
                    sound_velocity = 1.0e+12 # Default sound velocity value
                    KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(el.Properties.Id, sound_velocity))
                if sound_velocity <= 0.0:
                    raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, el.Properties.Id))
            # Get RESISTANCE
            if set_resistance:
                if el.Properties.Has(KratosCFD.RESISTANCE):
                    resistance = el.Properties.GetValue(KratosCFD.RESISTANCE)
                else:
                    resistance = 0.0 # Default resistance value = 0.0
                    KratosMultiphysics.Logger.PrintWarning('No \'RESISTANCE\' value found in Properties {0}. Setting default value {1}'.format(el.Properties.Id, resistance))
                if resistance < 0.0: # RESISTANCE = 0: fluid phase of the domain
                    raise Exception("RESISTANCE set to {0} in Properties {1}, positive number expected.".format(resistance, el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_density:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self._GetLocalMeshNodes())
        if set_viscosity:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self._GetLocalMeshNodes())
        if set_sound_velocity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self._GetLocalMeshNodes())
        if set_resistance:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCFD.RESISTANCE, resistance, self._GetLocalMeshNodes())

    def AddDofs(self):
        dofs_and_reactions_to_add = []
        if (not self.is_adjoint): # NS
            dofs_and_reactions_to_add.append("VELOCITY_X")
            dofs_and_reactions_to_add.append("VELOCITY_Y")
            dofs_and_reactions_to_add.append("VELOCITY_Z")
            dofs_and_reactions_to_add.append("PRESSURE")
            dofs_str = "NS"
        else: # ADJ
            dofs_and_reactions_to_add.append("VELOCITY_ADJ_X")
            dofs_and_reactions_to_add.append("VELOCITY_ADJ_Y")
            dofs_and_reactions_to_add.append("VELOCITY_ADJ_Z")
            dofs_and_reactions_to_add.append("PRESSURE_ADJ")
            dofs_str = "ADJ NS"
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid Topology Optimization " + dofs_str + " solver DOFs added correctly.")

    def SolveSolutionStep(self):
            problem_physiscs = self._GetTopologyOptimizationStage()
            if (problem_physiscs == 1):
                problem_phase_str = "FLUID"
            elif (problem_physiscs == 2):
                problem_phase_str = "ADJ-F"
            else: 
                problem_phase_str = "ERROR|"
            if self.IsNodesIdsGlobalToLocalDictionaryEmpty():
                raise RuntimeError("Executing 'SolveSolutionStep' of FluidTopologyOptimizationSolverMpi' with self.nodes_ids_global_to_local_partition_dictionary == { }")
            self.PrintPhysicsParametersUpdateStatus(problem_phase_str)
            self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " solution step...")
            # Call the base fluid solver to solve current time step
            is_converged = super().SolveSolutionStep()
            self.MpiPrint("--|" + problem_phase_str + "| ---> Step Solved!")
            return is_converged
    
    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.CloneTimeStep(new_time)
        if (not self.IsAdjoint()): # NS
            self.main_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP] += 1
        else: 
            self.main_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP] += 1
        return new_time
    
    def _CustomCloneTimeStep(self, new_time):
        if (not self.IsAdjoint()):
            self.main_model_part.CloneTimeStep(new_time)
        else:
            self.ShiftAdjointBuffer(adjoint_vars=[KratosMultiphysics.VELOCITY_ADJ, KratosMultiphysics.PRESSURE_ADJ])
            self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] = new_time

    def ShiftAdjointBuffer(self, adjoint_vars):
        for node in self._GetLocalMeshNodes():
            for var in adjoint_vars:
                for i in range(self.min_buffer_size-1, 0, -1):
                    node.SetSolutionStepValue(var, i, node.GetSolutionStepValue(var, i-1))

    def IsAdjoint(self):
        return self.is_adjoint
    
    def IsPhysics(self):
        return not self.IsAdjoint()
    
    def GetMainModelPart(self):
        return self.main_model_part
    
    def _SetTopologyOptimizationStage(self, top_opt_stage):
        self.GetComputingModelPart().ProcessInfo.SetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE, top_opt_stage)

    def _GetTopologyOptimizationStage(self):
        return self.GetComputingModelPart().ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
    
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
            source_mp = model_parts[0]
            modeler = KratosMultiphysics.ConnectivityPreserveModeler()
            element_name, condition_name = self.__GetElementAndConditionNames()
            modeler.GenerateModelPart(source_mp, self.main_model_part, element_name, condition_name)
            if _CheckIsDistributed():
                self.comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
                ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.main_model_part.GetRootModelPart(), self.comm)
                ParallelFillCommunicator.Execute()
                self.distributed_model_part_importer = physics_solver_distributed_model_part_importer

    def PrepareModelPart(self):
        # Call the base solver to do the PrepareModelPart
        # Note that his also calls the PrepareModelPart of the turbulence model
        super().PrepareModelPart()
        # Create the MPI communicators
        if _CheckIsDistributed():
            self.distributed_model_part_importer.CreateCommunicators()

    def _UpdateResistanceVariable(self, resistance):
        self.is_resistance_updated = True
        if isinstance(resistance, (int, float)): 
            for node in self._GetLocalMeshNodes():
             node.SetValue(KratosCFD.RESISTANCE, resistance)
        elif isinstance(resistance, (np.ndarray, list)): 
            for node in self._GetLocalMeshNodes():
                node.SetValue(KratosCFD.RESISTANCE, resistance[self.nodes_ids_global_to_local_partition_dictionary[node.Id]])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateResistanceVariable' : {type(resistance)}")
        
    def _UpdateConvectionVelocityVariable(self, convection_velocity):
        dim = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        for node in self._GetLocalMeshNodes():
            nodal_convection_velocity = convection_velocity[self.nodes_ids_global_to_local_partition_dictionary[node.Id]]
            node.SetValue(KratosCFD.CONVECTION_VELOCITY_X, nodal_convection_velocity[0])
            node.SetValue(KratosCFD.CONVECTION_VELOCITY_Y, nodal_convection_velocity[1])
            if (dim == 3):
                node.SetValue(KratosCFD.CONVECTION_VELOCITY_Z, nodal_convection_velocity[2])
    
    def InitializeSolutionStep(self):
        self.is_resistance_updated = False
        self._GetSolutionStrategy().InitializeSolutionStep()

    def PrintPhysicsParametersUpdateStatus(self, problem_phase_str):
        if (not self.is_resistance_updated):
            self.MpiPrint("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating RESISTANCE variable")

    def _CheckMaterialProperties(self):
        for node in self._GetLocalMeshNodes():
            self.MpiPrint("--|--> Resistance: " + node.GetValue(KratosCFD.RESISTANCE))
            break

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

    def SetNodesIdsGlobalToLocalDictionary(self, dict):
        self.nodes_ids_global_to_local_partition_dictionary = dict

    def IsNodesIdsGlobalToLocalDictionaryEmpty(self):
        return self.nodes_ids_global_to_local_partition_dictionary == {}
    
    def _GetLocalMeshNodes(self, mp = None):
        if mp is None:
            return self.GetMainModelPart().GetCommunicator().LocalMesh().Nodes
        else:
            return mp.GetCommunicator().LocalMesh().Nodes
        
    def __GetElementAndConditionNames(self):
        base_element_name = self.element_name
        base_condition_name = self.condition_name
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
        


