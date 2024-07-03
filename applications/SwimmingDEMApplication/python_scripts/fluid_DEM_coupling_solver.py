
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.SwimmingDEMApplication as KratosSDEM
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.SwimmingDEMApplication.CFD_DEM_coupling import ProjectionModule

def CreateSolver(model, custom_settings, projector):
    return FluidDEMSolver(model, custom_settings, projector)

class FluidDEMSolver(FluidSolver):

    ## FluidDEMSolver specific methods.
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step if we are not using a manufactured solution (step 0, read from input)
        if self.main_model_part.ProcessInfo[KratosSDEM.MANUFACTURED]:
            return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 2 >= self.GetMinimumBufferSize()
        else:
            return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()


    def _CreateScheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Cases in which the element manages the time integration
        if self.element_integrates_in_time:
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                domain_size,
                domain_size + 1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        # Cases in which a time scheme manages the time integration
        else:
            # Bossak time integration scheme
            if self.settings["time_scheme"].GetString() == "bossak":
                if self.settings["consider_periodic_conditions"].GetBool() == True:
                    scheme = KratosSDEM.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled(
                        self.settings["alpha"].GetDouble(),
                        domain_size,
                        KratosSDEM.PATCH_INDEX)
                else:
                    scheme = KratosSDEM.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled(
                        self.settings["alpha"].GetDouble(),
                        self.settings["move_mesh_strategy"].GetInt(),
                        domain_size)
            # BDF2 time integration scheme
            elif self.settings["time_scheme"].GetString() == "bdf2":
                scheme = KratosSDEM.BDF2TurbulentSchemeDEMCoupled()
            # Time scheme for steady state fluid solver
            elif self.settings["time_scheme"].GetString() == "steady":
                scheme = KratosSDEM.ResidualBasedSimpleSteadySchemeDEMCoupled(
                        self.settings["velocity_relaxation"].GetDouble(),
                        self.settings["pressure_relaxation"].GetDouble(),
                        domain_size)
            else:
                err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                err_msg += "Available options are: \"bossak\", \"bdf2\" and \"steady\""
                raise Exception(err_msg)

        return scheme

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            solution_strategy = self._CreateLinearStrategy()
        elif analysis_type == "non_linear":
            solution_strategy = self._CreateNewtonRaphsonStrategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return solution_strategy

    def SolveSolutionStep(self,projection_module):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep(projection_module)
        if not is_converged:
            msg  = "Fluid solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, msg)
        return is_converged

    def _CreateNewtonRaphsonStrategy(self):
        import KratosMultiphysics.SwimmingDEMApplication.residual_based_newton_raphson_strategy as NewtonRaphsonStrategy
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return NewtonRaphsonStrategy.ResidualBasedNewtonRaphsonStrategyPython(
            computing_model_part,
            time_scheme,
            convergence_criterion,
            builder_and_solver,
            self.settings["maximum_iterations"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())