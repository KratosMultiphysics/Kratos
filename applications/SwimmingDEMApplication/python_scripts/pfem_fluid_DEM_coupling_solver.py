import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_solver import PfemFluidSolver


def CreateSolver(model, parameters):
    return PfemFluidDEMcouplingSolver(model, parameters)

class PfemFluidDEMcouplingSolver(PfemFluidSolver):

    def __init__(self, model, parameters):

        super().__init__(model, parameters)

    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("SwimmingDEM", self.main_model_part.SetBufferSize(self.settings["buffer_size"].GetInt()))
        
        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        self.fluid_solver = KratosPfemFluid.TwoStepVPDEMcouplingStrategy(self.computing_model_part,
                                                              self.velocity_linear_solver,
                                                              self.pressure_linear_solver,
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              self.settings["velocity_tolerance"].GetDouble(),
                                                              self.settings["pressure_tolerance"].GetDouble(),
                                                              self.settings["maximum_pressure_iterations"].GetInt(),
                                                              self.settings["time_order"].GetInt(),
                                                              self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION])

        echo_level = self.settings["echo_level"].GetInt()

        # Set echo_level
        self.fluid_solver.SetEchoLevel(echo_level)

        # Check if everything is assigned correctly
        self.fluid_solver.Check()

    def AddVariables(self):

        super().AddVariables()

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_OLD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_RATE)
