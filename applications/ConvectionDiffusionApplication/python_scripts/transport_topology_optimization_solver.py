# Importing the Kratos Library
import KratosMultiphysics

import numpy as np #import the numpy library
import scipy as sp #import the scipy library


# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Importing the base class
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_transient_solver import ConvectionDiffusionTransientSolver
import KratosMultiphysics.FluidDynamicsApplication

def CreateSolver(model, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return TransportTopologyOptimizationSolver(model, solver_settings, is_adjoint_solver=isAdjointSolver)

class TransportTopologyOptimizationSolver(ConvectionDiffusionTransientSolver):

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
        super().__init__(model,custom_settings)
        self._DefineAdjointSolver(is_adjoint_solver)
        self._DefineElementsAndConditions()
        # self._InitializePhysicsParameters()
        print_str = "Construction of TransportTopologyOptimizationSolver "
        if self.IsAdjoint():
            print_str += "for Adjoint problem "
        print_str +=  "finished."
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, print_str)

    def _DefineAdjointSolver(self, isAdjointSolver):
        self.is_adjoint = isAdjointSolver

    def _DefineElementsAndConditions(self):
        self.element_name = "TransportTopologyOptimizationElement"
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
        num_nodes_elements = len(self.main_model_part.Elements[1].GetNodes())
        name_string = f"{self.element_name}{domain_size}D{num_nodes_elements}N"
        self.settings["element_replace_settings"]["element_name"].SetString(name_string)
        ## Conditions
        num_nodes_conditions = domain_size
        name_string = f"{self.condition_name}{domain_size}D{num_nodes_conditions}N"
        self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]

    def AddVariables(self): 
        # Add Transport Topology Optimization Variables 
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
        print("--|--> Conductivity:", self.main_model_part.Nodes[1].GetValue(KratosMultiphysics.CONDUCTIVITY))
        print("--|--> Decay:", self.main_model_part.Nodes[1].GetValue(KratosCD.DECAY))
        print("--|--> Convection Coefficient:", self.main_model_part.Nodes[1].GetValue(KratosMultiphysics.CONVECTION_COEFFICIENT))  
    
    def _SetTimeSchemeBufferSize(self):
        self.settings["time_scheme"].SetString("bdf2")
        self.min_buffer_size = 1

    def PrepareModelPart(self):
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.CONDUCTIVITY, 0.0, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCD.DECAY, 0.0, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.CONVECTION_COEFFICIENT, 0.0, self.main_model_part.Nodes)
        if not self.is_restarted():
            # Import material properties
            materials_imported = self.import_materials()
            if materials_imported:
                KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were successfully imported.")
            else:
                KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]:: ", "Materials were not imported.")

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()
            self._set_and_fill_buffer()

        if (self.settings["echo_level"].GetInt() > 0):
            KratosMultiphysics.Logger.PrintInfo(self.model)
        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionSolver]::", "ModelPart prepared for Solver.")

    def _UpdateConductivityVariable(self, conductivity):
        self.is_conductivity_updated = True
        mp = self.GetComputingModelPart()
        if isinstance(conductivity, (int, float)): 
            for node in mp.Nodes:
             node.SetValue(KratosMultiphysics.CONDUCTIVITY, conductivity)
        elif isinstance(conductivity, (np.ndarray, list)): 
            for node in mp.Nodes:
                node.SetValue(KratosMultiphysics.CONDUCTIVITY, conductivity[node.Id-1])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateConductivityVariable' : {type(conductivity)}")
        
    def _UpdateDecayVariable(self, decay):
        self.is_decay_updated = True
        mp = self.GetComputingModelPart()
        if isinstance(decay, (int, float)): 
            for node in mp.Nodes:
             node.SetValue(KratosCD.DECAY, decay)
        elif isinstance(decay, (np.ndarray, list)): 
            for node in mp.Nodes:
                node.SetValue(KratosCD.DECAY, decay[node.Id-1])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateDecayVariable' : {type(decay)}")

    def _UpdateConvectionCoefficientVariable(self, convection_coefficient):
        self.is_convection_coefficient_updated = True
        mp = self.GetComputingModelPart()
        if isinstance(convection_coefficient, (int, float)): 
            for node in mp.Nodes:
             node.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, convection_coefficient)
        elif isinstance(convection_coefficient, (np.ndarray, list)): 
            for node in mp.Nodes:
                node.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, convection_coefficient[node.Id-1])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateResistanceVariable' : {type(convection_coefficient)}")
        
    def _UpdateTransportSourceVariable(self, transport_source):
        self.is_transport_source_updated = True
        mp = self.GetComputingModelPart()
        if isinstance(transport_source, (int, float)): 
            for node in mp.Nodes:
             node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, transport_source)
        elif isinstance(transport_source, (np.ndarray, list)): 
            for node in mp.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, transport_source[node.Id-1])
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
        mp = self.GetComputingModelPart()
        for node in mp.Nodes:
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
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = 1.0
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 0.0
        transport_top_opt_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return transport_top_opt_scheme
    
    def SolveSolutionStep(self):
        problem_physiscs = self._GetTopologyOptimizationStage()
        if (problem_physiscs == 1):
            problem_phase_str = "TRANS"
        elif (problem_physiscs == 2):
            problem_phase_str = "ADJ-T"
        else: 
            problem_phase_str = "ERROR|"
        self.PrintPhysicsParametersUpdateStatus(problem_phase_str)
        print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " solution step...")
        # Call the base fluid solver to solve current time step
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        print("--|" + problem_phase_str + "| ---> Step Solved!")
        return is_converged
    
    def AdvanceInTime(self, current_time):
        if (not self.IsAdjoint()): # T
            self.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_T_STEP] += 1
            new_time =  super().AdvanceInTime(current_time)
        else: #ADJ
            print("[WARNING] The Adjoint problem should go backward buth this has not yet been implemented. Since for now it is steady, we advance in time even if it is unnecessary.")
            dt = self.ComputeDeltaTime()
            # new_time = current_time - dt
            new_time = current_time + dt
            self.main_model_part.CloneTimeStep(new_time)
            # print("\nASK HOW TO HANDLE THIS!!!\n")
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
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
    
    def InitializeSolutionStep(self):
        # if (self.IsAdjoint()):
        #     print("strat", self._GetSolutionStrategy())
        #     print("f_w:", self.GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS))
        #     input("adjoint solver initialize sol step")
        self._GetSolutionStrategy().InitializeSolutionStep()

    def ImportModelPart(self, model_parts=None):
        if (self.IsPhysics()):
            # Call the fluid solver to import the model part from the mdpa
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])
        else:
            transport_mp = model_parts[0]
            self.main_model_part = transport_mp
    
    def InitializeSolutionStep(self):
        self.is_conductivity_updated = False
        self.is_decay_updated = False
        self.is_convection_coefficient_updated = False
        self.is_transport_source_updated = True
        self._GetSolutionStrategy().InitializeSolutionStep()

    def PrintPhysicsParametersUpdateStatus(self, problem_phase_str):
        if (not self.is_conductivity_updated):
            print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating CONDUCTIVITY variable")
        if (not self.is_decay_updated):
            print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating DECAY variable")
        if (not self.is_convection_coefficient_updated):
            print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating CONVECTION_COEFFICIENT variable")
        if (not self.is_transport_source_updated):
            print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating TRANSPORT SOURCE variable")
