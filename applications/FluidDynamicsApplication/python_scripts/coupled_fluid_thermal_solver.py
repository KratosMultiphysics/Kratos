from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff

def CreateSolver(main_model_part, custom_settings):
    return CoupledFluidThermalSolver(main_model_part, custom_settings)

class CoupledFluidThermalSolver(object):

    def __init__(self, main_model_part, custom_settings):
        
        self.main_model_part = main_model_part
                
        default_settings = KratosMultiphysics.Parameters("""
            {
                "solver_type" : "ThermallyCoupled",
                "fluid_solver_settings": {
                        "solver_type": "navier_stokes_solver_vmsmonolithic",
                        "model_import_settings": {
                                "input_type": "mdpa",
                                "input_filename": "unknown_name"
                        }
                },
                "thermal_solver_settings": {
                        "solver_type": "Transient",
                        "analysis_type": "linear",
                        "model_import_settings": {
                                "input_type": "use_input_model_part",
                                "input_filename": "unknown_name"
                        },
                        "computing_model_part_name": "Thermal",
                        "material_import_settings": {
                                "materials_filename": "ThermicMaterials.json"
                        },
                        "convection_diffusion_variables": {
                                "density_variable": "DENSITY",
                                "diffusion_variable": "CONDUCTIVITY",
                                "unknown_variable": "TEMPERATURE",
                                "volume_source_variable": "HEAT_FLUX",
                                "surface_source_variable": "FACE_HEAT_FLUX",
                                "projection_variable": "PROJECTED_SCALAR1",
                                "convection_variable": "CONVECTION_VELOCITY",
                                "mesh_velocity_variable": "MESH_VELOCITY",
                                "transfer_coefficient_variable": "",
                                "velocity_variable": "VELOCITY",
                                "specific_heat_variable": "SPECIFIC_HEAT",
                                "reaction_variable": "REACTION_FLUX"
                        }
                }
            }
            """)

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.main_model_part, self.settings["fluid_solver_settings"],"OpenMP")
        
        
        self.thermal_model_part =  KratosMultiphysics.ModelPart("thermal_model_part")
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(self.main_model_part, self.thermal_model_part, "Element2D3N", "Condition2D2N")
        
        import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.thermal_model_part,self.settings["thermal_solver_settings"],"OpenMP")
        

    def AddVariables(self):
        self.fluid_solver.AddVariables()
        
        #this is a HACK: TODO: cleanup so that this is not needed in the future 
        #the problem is that variables are not added to the fluid_solver.MainModelPart who is in charge of creating the nodes
        self.tmp = self.thermal_solver.main_model_part
        self.thermal_solver.main_model_part = self.fluid_solver.main_model_part
        
        self.thermal_solver.AddVariables()

        #this is a HACK: TODO: cleanup so that this is not needed in the future 
        self.thermal_solver.main_model_part = self.tmp
        

    def ImportModelPart(self):
        self.fluid_solver.ImportModelPart()
        
        #here cloning the fluid modelpart to thermal_model_part so that the nodes are shared
        #convection_diffusion_settings = self.thermal_solver.GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(self.main_model_part, self.thermal_model_part, "Element2D3N", "Condition2D2N")
        #self.thermal_solver.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)


        self.thermal_solver.ImportModelPart()
        
    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        return self.fluid_solver.ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()

    def SaveRestart(self):
        pass #one should write the restart file here

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt

        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        self.thermal_solver.SolveSolutionStep()
        
    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()


