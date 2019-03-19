from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
#import KratosMultiphysics.SolidMechanicsApplication as KratosSolidMechanics

import pfem_fluid_solver as BaseSolver

def CreateSolver(model, parameters):
    return PfemFluidExplicitSolver(model, parameters)

class PfemFluidExplicitSolver(BaseSolver.PfemFluidSolver):

    def __init__(self, main_model_part, custom_settings):

         # Construct the base solver.
        super(PfemFluidExplicitSolver,self).__init__(model,parameters)

        print("Construction of Pfem Fluid Explicit Solver finished.")


    def Initialize(self):

        print("::[Pfem Fluid Explicit Solver]:: -START-")

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        #mechanical_scheme = KratosSolidMechanics.ExplicitCentralDifferencesScheme(self.settings["max_delta_time"].GetDouble(),
        #                                                                          self.settings["fraction_delta_time"].GetDouble(),
        #                                                                          self.settings["time_step_prediction_level"].GetDouble(),
        #                                                                          self.settings["rayleigh_damping"].GetBool())

        mechanical_scheme = KratosPfemFluid.FirstOrderForwardEulerScheme(1.0e-4,
                                                                         1.0,
                                                                         0,
                                                                         0)
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])
        #self.fluid_solver = KratosPfemFluid.ExplicitStrategy(self.computing_model_part,
        #                                                     mechanical_scheme,
        #                                                     linear_solver,
        #                                                     self.settings["compute_reactions"].GetBool(),
        #                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
        #                                                     self.settings["move_mesh_flag"].GetBool())

        self.fluid_solver = KratosPfemFluid.ExplicitStrategy(self.computing_model_part,
                                                             mechanical_scheme,
                                                             linear_solver,
                                                             False,
                                                             True,
                                                             True)
        # Set echo_level
        self.fluid_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set initialize flag
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.mechanical_solver.SetInitializePerformedFlag(True)

        # Check if everything is assigned correctly
        self.fluid_solver.Check()


        print("::[Pfem Fluid  Explicit Solver]:: -END- ")
