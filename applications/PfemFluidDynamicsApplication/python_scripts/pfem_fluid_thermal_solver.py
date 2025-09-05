import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
from KratosMultiphysics.python_solver import PythonSolver

from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_fluid_solver as BaseSolver

def CreateSolver(model, parameters):
    return PfemFluidThermalSolver(model, parameters)

class PfemFluidThermalSolver(BaseSolver.PfemFluidSolver):

    def Initialize(self):

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        self.fluid_solver = KratosPfemFluid.TwoStepVPThermalStrategy(self.computing_model_part,
                                                               self.velocity_linear_solver,
                                                               self.pressure_linear_solver,
                                                               self.settings["reform_dofs_at_each_step"].GetBool(),
                                                               self.settings["velocity_tolerance"].GetDouble(),
                                                               self.settings["pressure_tolerance"].GetDouble(),
                                                               self.settings["maximum_pressure_iterations"].GetInt(),
                                                               self.settings["time_order"].GetInt(),
                                                               self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION])

        # Set echo_level
        echo_level = self.settings["echo_level"].GetInt()
        self.fluid_solver.SetEchoLevel(echo_level)

        # Self initialize strategy
        self.fluid_solver.Initialize()

        # Check if everything is assigned correctly
        self.fluid_solver.Check() #TODO: This must be done in the Check function