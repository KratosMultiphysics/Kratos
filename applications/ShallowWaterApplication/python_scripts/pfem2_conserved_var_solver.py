from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.pfem2_primitive_var_solver import Pfem2PrimitiveVarSolver

def CreateSolver(model, custom_settings):
    return Pfem2ConservedVarSolver(model, custom_settings)

class Pfem2ConservedVarSolver(Pfem2PrimitiveVarSolver):
    def __init__(self, model, custom_settings):
        super(Pfem2ConservedVarSolver,self).__init__(model, custom_settings)

        # Set the element and condition names for the replace settings
        self.element_name = "ConservedVarElement"
        self.condition_name = "Condition"
        self.min_buffer_size = 2

        # Pfem2 settings
        self.settings["pfem2_settings"]["convection_vector_variable"].SetString("MOMENTUM")

    def AddVariables(self):
        super(Pfem2ConservedVarSolver,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

        KM.Logger.PrintInfo("::[Pfem2ConservedVarSolver]::", "Shallow water solver DOFs added correctly.")

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Solve equations on mesh
            is_converged = self.solver.SolveSolutionStep()
            # Compute free surface
            SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.main_model_part)
            # Compute velocity
            SW.ShallowWaterUtilities().ComputeVelocity(self.main_model_part)
            # Print particles if needed
            if self.print_particles:
                self.lagrangian_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
                self.lagrangian_model_part.ProcessInfo[KratosMultiphysics.TIME] = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
                self.moveparticles.ExecuteParticlesPrintingTool(self.lagrangian_model_part, self.filter_factor)

            return is_converged
