# Importing the Kratos Library
import KratosMultiphysics
from importlib import import_module

import numpy as np #import the numpy library
import scipy as sp #import the scipy library


# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import NavierStokesMonolithicSolver

def CreateSolver(model, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return FluidTopologyOptimizationSolver(model, solver_settings, is_adjoint_solver=isAdjointSolver)

class FluidTopologyOptimizationSolver(NavierStokesMonolithicSolver):

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
        super().__init__(model,custom_settings)
        self._DefineAdjointSolver(is_adjoint_solver)
        self._DefineElementsAndConditions()
        print_str = "Construction of FluidTopologyOptimizationSolver "
        if self.IsAdjoint():
            print_str += "for Adjoint problem "
        print_str +=  "finished."
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, print_str)

    def _DefineAdjointSolver(self, isAdjointSolver):
        self.is_adjoint = isAdjointSolver

    def _DefineElementsAndConditions(self):
        if (self.element_name != "FluidTopologyOptimizationElement"):
            print("[WARNING]", self.__class__.__name__, "element_name: \'", self.element_name, "\' is not compatible with FluidTopologyOptimization. Its value has been reset to default value: \' FluidTopologyOptimizationElement \'")
            self.element_name = "FluidTopologyOptimizationElement"
        self.condition_name = "NavierStokesWallCondition"
        self.element_integrates_in_time = True
        
    def _SetFormulation(self):
        super()._SetFormulation()
        self.element_has_nodal_properties = True
        self.non_historical_nodal_properties_variables_list.append(KratosCFD.RESISTANCE)
        self.non_historical_nodal_properties_variables_list.append(KratosMultiphysics.CONVECTION_COEFFICIENT)

    def AddVariables(self):  
        #Add parent class variables
        super().AddVariables()    
        # Add Adjoint-NS Variables
        # self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ACCELERATION_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.MESH_VELOCITY_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.BODY_FORCE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        # eventual transport coupling
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_ADJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_ADJ_GRADIENT)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid Topology Optimization ADJ-NS solver variables added correctly.")          

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bossak":
            print("[WARNING]", self.__class__.__name__, "time scheme_type: \'", scheme_type, "\' is not compatible with the current implementation of FluidTopologyOptimization. Its value has been reset to default value: \' steady \'")
            self.settings["time_scheme"].SetString("steady")
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        elif scheme_type == "bdf2":
            print("[WARNING]", self.__class__.__name__, "time scheme_type: \'", scheme_type, "\' is not compatible with the current implementation of FluidTopologyOptimization. Its value has been reset to default value: \' steady \'")
            self.settings["time_scheme"].SetString("steady")
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        elif scheme_type == "steady":
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"bossak\", \"bdf2\" or \"steady\".\n"
            raise Exception(msg)
        # TO DO: should be set to steady, but due to its implementation works only with 
        # scheme:bdf2 
        # min buffer size = 1
        self.settings["time_scheme"].SetString("bdf2")
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
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        if set_viscosity:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
        if set_sound_velocity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)
        if set_resistance:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCFD.RESISTANCE, resistance, self.main_model_part.Nodes)

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
        self.PrintPhysicsParametersUpdateStatus(problem_phase_str)
        print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " solution step...")
        # Call the base fluid solver to solve current time step
        is_converged = super().SolveSolutionStep()
        print("--|" + problem_phase_str + "| ---> Step Solved!")
        return is_converged
    
    def AdvanceInTime(self, current_time):
        if (not self.IsAdjoint()): # NS
            self.main_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP] += 1
            new_time =  super().AdvanceInTime(current_time)
        else: #ADJ
            print("[WARNING] The Adjoint problem should go backward buth this has not yet been implemented. Since for now it is steady, we advance in time even if it is unnecessary.")
            dt = self._ComputeDeltaTime()
            # new_time = current_time - dt
            new_time = current_time + dt
            self.main_model_part.CloneTimeStep(new_time)
            # print("\nASK HOW TO HANDLE THIS!!!\n")
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            self.main_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP] += 1
        return new_time
    
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
    
    def ImportModelPart(self, model_parts=None):
        if (self.IsPhysics()):
            # Call the fluid solver to import the model part from the mdpa
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])
        else:
            fluid_mp = model_parts[0]
            self.main_model_part = fluid_mp

    def _UpdateResistanceVariable(self, resistance):
        self.is_resistance_updated = True
        mp = self.GetComputingModelPart()
        if isinstance(resistance, (int, float)): 
            for node in mp.Nodes:
             node.SetValue(KratosCFD.RESISTANCE, resistance)
        elif isinstance(resistance, (np.ndarray, list)): 
            for node in mp.Nodes:
                node.SetValue(KratosCFD.RESISTANCE, resistance[node.Id-1])
        else:
            raise TypeError(f"Unsupported input type in '_UpdateResistanceVariable' : {type(resistance)}")
    
    def InitializeSolutionStep(self):
        self.is_resistance_updated = False
        self._GetSolutionStrategy().InitializeSolutionStep()

    def PrintPhysicsParametersUpdateStatus(self, problem_phase_str):
        if (not self.is_resistance_updated):
            print("--|" + problem_phase_str + "| ---> Top. Opt. solution: Solve " + problem_phase_str + " without updating RESISTANCE variable")

    def _CheckMaterialProperties(self):
        print("--|--> Resistance:", self.main_model_part.Nodes[1].GetValue(KratosCFD.RESISTANCE))

        


