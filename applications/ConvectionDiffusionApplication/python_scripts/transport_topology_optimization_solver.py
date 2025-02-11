# Importing the Kratos Library
import KratosMultiphysics


# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Importing the base class
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_transient_solver import ConvectionDiffusionTransientSolver
import KratosMultiphysics.FluidDynamicsApplication

def CreateSolver(model, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return TransportTopologyOptimizationSolver(model, solver_settings, isAdjointSolver)

class TransportTopologyOptimizationSolver(ConvectionDiffusionTransientSolver):
    def __init__(self, model, custom_settings, isAdjointSolver = False):
        super().__init__(model,custom_settings)
        self._DefineAdjointSolver(isAdjointSolver)
        self._DefineElementsAndConditions(custom_settings)
        self._InitializePhysicsProperties()
        print_str = "Construction of TransportTopologyOptimizationSolver "
        if isAdjointSolver:
            print_str += "for Adjoint problem "
        print_str +=  "finished."
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, print_str)

    def _DefineAdjointSolver(self, isAdjointSolver):
        self.is_adjoint = isAdjointSolver

    def _DefineConvectionVelocity(self, const_velocity, vel=[0,0,0]):
        self.is_constant_velocity = const_velocity
        self.velocity = vel

    def _InitializePhysicsProperties(self):
        self.decay = 0
        self.is_constant_velocity = False
        self.velocity = [0,0,0]
    
    def _DefineNodalProperties(self, decay):
        self.decay = decay

    def _DefineElementsAndConditions(self,settings):
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
        ## Note that we check for the elements that require substitution to allow for custom elements
        domain_size = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        num_nodes_elements = domain_size + 1
        name_string = f"{self.element_name}{domain_size}D{num_nodes_elements}N"
        self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        num_nodes_conditions = domain_size
        name_string = f"{self.condition_name}{domain_size}D{num_nodes_conditions}N"
        self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]

    def AddVariables(self): 
        # Add Fluid Transport Topology Optimization Variables 
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.DECAY)
        # if (not self.IsAdjoint()):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        # else:
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.TEMPERATURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.HEAT_FLUX_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCD.FACE_HEAT_FLUX_ADJ)
        self._DefineSettings()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid Transport Topology Optimization T and T_ADJ solver variables added correctly.")
        
    def _DefineSettings(self):
        super().AddVariables()
    
    def _CheckMaterialProperties(self):
        print("Conductivity:", self.main_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY))
        print("Decay:", self.main_model_part.Nodes[1].GetSolutionStepValue(KratosCD.DECAY))

    def _CheckProperties(self):
        self._CheckMaterialProperties()
    
    def _SetTimeSchemeBufferSize(self):
        self.settings["time_scheme"].SetString("bdf2")
        self.min_buffer_size = 1
        # self._SetUpSteadySimulation()

    def PrepareModelPart(self):
        self._SetDecay()
        self._SetConvectiveVelocity()
        super().PrepareModelPart()
        # self._CheckProperties()    
        
    def _SetDecay(self):
        mp = self.GetComputingModelPart()
        for el in mp.Elements:
            el.Properties.SetValue(KratosCD.DECAY, self.decay)

    def _SetConvectiveVelocity(self):
        if (self.IsConstantVelocity()):
            self._SetConstantConvectiveVelocity()
        else:
            self._SetNonConstantConvectiveVelocity()
    
    def _SetConstantConvectiveVelocity(self):
        mp = self.GetComputingModelPart()
        for node in mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, KratosMultiphysics.Vector(self.velocity))

    def _SetNonConstantConvectiveVelocity(self):
        pass

    def IsConstantVelocity(self):
        return self.is_constant_velocity
    
    def AddDofs(self):
        dofs_and_reactions_to_add = []
        if (not self.IsPhysics()): # T
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
        

