import KratosMultiphysics
import KratosMultiphysics.SPHApplication as SPHApplication
from KratosMultiphysics.SPHApplication.sph_solver import SPHSolver

def CreateSolver(model, custom_settings):
    return SPHTransientSolver(model, custom_settings)


class SPHTransientSolver(SPHSolver):

    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 2

        theta = self.settings["transient_parameters"]["theta"].GetDouble()
        if not 0.0 < theta <= 1.0:
            raise ValueError('"transient_parameters.theta" must be in the interval (0, 1]. ' f"Got {theta}.")
        if not self.settings["strain_dofs"].GetBool():
            raise ValueError('The transient SPH solver requires the mixed formulation with "strain_dofs" set to true.')

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
            "time_integration_method" : "implicit",
            "strain_dofs" : true,
            "transient_parameters" : {
                "theta"    : 0.5
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super().AddVariables()
        self._add_dynamic_variables()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Variables ADDED")

    def _CreateScheme(self):
        # Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = self.settings["transient_parameters"]["theta"].GetDouble()

        # Time discretization is assembled by the element; the scheme only updates the DOFs.
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        sph_scheme = self._GetScheme()
        sph_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        strategy = SPHApplication.MixedSphResidualBasedNewtonRaphsonStrategy(
            computing_model_part,
            sph_scheme,
            sph_convergence_criterion,
            builder_and_solver,
            self.settings["max_iteration"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())
        strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings["use_old_stiffness_in_first_iteration"].GetBool())
        return strategy