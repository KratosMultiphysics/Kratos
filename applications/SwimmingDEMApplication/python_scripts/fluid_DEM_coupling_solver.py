from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.SwimmingDEMApplication as KratosSDEM
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return FluidDEMSolver(model, custom_settings)

class FluidDEMSolver(FluidSolver):

    ## FluidDEMSolver specific methods.

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
                scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                        self.settings["velocity_relaxation"].GetDouble(),
                        self.settings["pressure_relaxation"].GetDouble(),
                        domain_size)
            else:
                err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                err_msg += "Available options are: \"bossak\", \"bdf2\" and \"steady\""
                raise Exception(err_msg)

        return scheme