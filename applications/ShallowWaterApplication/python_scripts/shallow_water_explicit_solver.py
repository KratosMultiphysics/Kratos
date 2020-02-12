from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library and applications
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return ShallowWaterExplicitSolver(model, custom_settings)

class ShallowWaterExplicitSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(ShallowWaterExplicitSolver, self).__init__(model, settings)
        self.element_name = "ConservedElement"
        self.condition_name = "Condition"
    
    def AddVariables(self):
        super(ShallowWaterExplicitSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(SW.MOMENTUM_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(SW.MOMENTUM_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(SW.HEIGHT_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(SW.HEIGHT_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_MASS)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)
        
    def Initialize(self):
        # The time step utility needs the NODAL_H
        KM.FindNodalHProcess(self.GetComputingModelPart()).Execute()
        self.EstimateDeltaTimeUtility = SW.EstimateDtShallow(self.GetComputingModelPart(), self.settings["time_stepping"])

        self.solver = SW.RungeKuttaStrategy(
            self.GetComputingModelPart(),
            self.GetComputingModelPart().ProcessInfo[KM.DOMAIN_SIZE],
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())
        
        self.main_model_part.ProcessInfo.SetValue(KM.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        self.solver.SetEchoLevel(max(0, self.echo_level-1))
        self.solver.Check()

        self.solver.Initialize()

        KM.Logger.PrintInfo("::[ShallowWaterExplicitSolver]::", "Mesh stage solver initialization finished")