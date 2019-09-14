from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_fluid_solver as BaseSolver

def CreateSolver(model, parameters):
    return PfemFluidExplicitSolver(model, parameters)

class PfemFluidExplicitSolver(BaseSolver.PfemFluidSolver):
    """The base class for the PfemFluidExplicitSolver-classes
    """
    def __init__(self, main_model_part, custom_settings):
        """The constructor of the PfemFluidExplicitSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        main_model_part -- The ModelPart to be used
        custom_settings -- The settings for the solver
        """
         # Construct the base solver.
        super(PfemFluidExplicitSolver,self).__init__(model,parameters)
        self.KratosPrintInfo("Construction of Pfem Fluid Explicit Solver finished.")


    def Initialize(self):
        """This function initializes the PfemFluidExplicitSolver
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self.KratosPrintInfo("::[Pfem Fluid Explicit Solver]:: -START-")

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = KratosPfemFluid.FirstOrderForwardEulerScheme(1.0e-4, 1.0, 0, 0)
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])
        self.fluid_solver = KratosPfemFluid.ExplicitStrategy(self.computing_model_part,
                                                             mechanical_scheme,
                                                             linear_solver,
                                                             False,
                                                             True,
                                                             True)
        # Set echo_level
        self.fluid_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set initialize flag
        if self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True:
            self.mechanical_solver.SetInitializePerformedFlag(True)

        # Check if everything is assigned correctly
        self.fluid_solver.Check()
        self.KratosPrintInfo("::[Pfem Fluid  Explicit Solver]:: -END- ")


    def KratosPrintInfo(self, message):
        """This function prints info on screen
        """
        KratosMultiphysics.Logger.Print(message, label="")
        KratosMultiphysics.Logger.Flush()