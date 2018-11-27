from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication", "MeshMovingApplication")
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# other imports
from ale_fluid_solver import AleFluidSolver
import python_solvers_wrapper_fluid


def CreateSolver(model, solver_settings, parallelism):
    return NavierStokesAleFluidSolver(model, solver_settings, parallelism)


class NavierStokesAleFluidSolver(AleFluidSolver):
    def _CreateFluidSolver(self, solver_settings, parallelism):
        return python_solvers_wrapper_fluid.CreateSolverByParameters(
            self.model, solver_settings, parallelism)

    def _SelectMeshVelocityComputationSettings(self):
        '''Selecting the time-integration for the MESH_VELOCITY to be consistent
        with the fluid time-integration.
        By now the parameters of the FluidSolver have been validated, which means
        that the time-integration method used by the fluid can be queried
        '''
        mesh_vel_comp_settings = self.settings["mesh_velocity_computation"]
        fluid_settings = self.settings["fluid_solver_settings"]

        fluid_solver_type = fluid_settings["solver_type"].GetString()
        if fluid_solver_type == "monolithic" or fluid_solver_type == "Monolithic":
            if fluid_settings.Has("time_scheme"):
                time_scheme_fluid = fluid_settings["time_scheme"].GetString()
                alpha_fluid = fluid_settings["alpha"].GetDouble()
                if mesh_vel_comp_settings.Has("time_scheme"):
                    time_scheme_mesh_vel = mesh_vel_comp_settings["time_scheme"].GetString()
                    if time_scheme_fluid != time_scheme_mesh_vel and self.is_printing_rank:
                        info_msg  = '"time_scheme" of the fluid (' + time_scheme_fluid
                        info_msg += ') is different from the\n"time_scheme" used for the '
                        info_msg += 'computation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                else:
                    mesh_vel_comp_settings.AddValue("time_scheme", fluid_settings["time_scheme"])
                    if self.is_printing_rank:
                        info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncomputation of the '
                        info_msg += 'mesh-velocity to "' + time_scheme_fluid + '" to be consistent with the\n'
                        info_msg += '"time_scheme" of the fluid'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                if not mesh_vel_comp_settings["time_scheme"].GetString().startswith("bdf"):
                    if mesh_vel_comp_settings.Has("alpha_m"):
                        alpha_mesh_vel = mesh_vel_comp_settings["alpha_m"].GetDouble()
                        if abs(alpha_fluid-alpha_mesh_vel) > 1e-12 and self.is_printing_rank:
                            info_msg  = '"alpha" of the fluid (' + str(alpha_fluid)
                            info_msg += ') is different from the\n"alpha_m" used for the '
                            info_msg += 'computation of the mesh-velocity (' + str(alpha_mesh_vel) + ')'
                            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                    else:
                        mesh_vel_comp_settings.AddValue("alpha_m", fluid_settings["alpha"])
                        if self.is_printing_rank:
                            info_msg  = 'setting "alpha_m" of the mesh-solver for the\ncomputation of the '
                            info_msg += 'mesh-velocity to "' + str(alpha_fluid) + '" to be consistent\nwith '
                            info_msg += '"alpha" of the fluid'
                            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        elif fluid_solver_type == "fractional_step" or fluid_solver_type == "FractionalStep":
            # currently fractional step always uses BDF2
            if mesh_vel_comp_settings.Has("time_scheme"):
                time_scheme_mesh_vel = mesh_vel_comp_settings["time_scheme"].GetString()
                if time_scheme_mesh_vel != "bdf2" and self.is_printing_rank:
                    info_msg  = '"time_scheme" of the fluid (bdf2) '
                    info_msg += 'is different from the\n"time_scheme" used for the '
                    info_msg += 'computation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                    KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
            else:
                mesh_vel_comp_settings.AddEmptyValue("time_scheme").SetString("bdf2")
                if self.is_printing_rank:
                    info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncomputation of the '
                    info_msg += 'mesh-velocity to "bdf2" to be consistent with the\n'
                    info_msg += '"time_scheme" of the fluid'
                    KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        if not mesh_vel_comp_settings.Has("time_scheme"):
            mesh_vel_comp_settings.AddEmptyValue("time_scheme").SetString("UNSPECIFIED")
            warn_msg  = 'unknown "solver_type" of the fluid-solver, therefore '
            warn_msg += 'no automatic selection\nof "time_scheme" for the computation '
            warn_msg += 'of the mesh-velocity performed'
            KratosMultiphysics.Logger.PrintWarning("::[ALEFluidSolver]::", warn_msg)
