from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# other imports
import KratosMultiphysics.FluidDynamicsApplication.python_solvers_wrapper_fluid as fluid_solvers_wrapper
from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

have_mesh_moving = CheckIfApplicationsAvailable("MeshMovingApplication")
if have_mesh_moving:
    from KratosMultiphysics.MeshMovingApplication.ale_fluid_solver import AleFluidSolver
else:
    raise Exception("In importing the NavierStokesAleFluidSolver: The solver requires the MeshMovingApplication, but this application is not available.")


def CreateSolver(model, solver_settings, parallelism):
    return NavierStokesAleFluidSolver(model, solver_settings, parallelism)


class NavierStokesAleFluidSolver(AleFluidSolver):
    def _CreateFluidSolver(self, solver_settings, parallelism):
        return fluid_solvers_wrapper.CreateSolverByParameters(
            self.model, solver_settings, parallelism)

    @classmethod
    def _ManipulateFluidSolverSettingsForReactionsComputation(cls, fluid_solver_settings):
        if fluid_solver_settings.Has("compute_reactions"):
            if not fluid_solver_settings["compute_reactions"].GetBool():
                fluid_solver_settings["compute_reactions"].SetBool(True)
                warn_msg  = '"compute_reactions" is switched off for the fluid-solver, switching it on!'
                KM.Logger.PrintWarning("::[NavierStokesAleFluidSolver]::", warn_msg)
        else:
            fluid_solver_settings.AddEmptyValue("compute_reactions").SetBool(True)
            info_msg = 'Setting "compute_reactions" to true for the fluid-solver'
            KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)

    @classmethod
    def _ManipulateMeshMotionSolverSettingsForMeshVelocityComputation(cls, fluid_solver_settings, mesh_motion_solver_settings):
        # 1. Ensure that the MESH_VELOCITY is computed
        if mesh_motion_solver_settings.Has("calculate_mesh_velocity"):
            if not mesh_motion_solver_settings["calculate_mesh_velocity"].GetBool():
                mesh_motion_solver_settings["calculate_mesh_velocity"].SetBool(True)
                KM.Logger.PrintWarning("::[NavierStokesAleFluidSolver]::", 'Mesh velocity calculation was deactivated. Switching "calculate_mesh_velocity" on')
        else:
            mesh_motion_solver_settings.AddEmptyValue("calculate_mesh_velocity").SetBool(True)

        # 2. Selecting the time-integration for the MESH_VELOCITY to be consistent with the fluid time-integration.
        #    By now the parameters of the FluidSolver have been validated, which means
        #    that the time-integration method used by the fluid can be queried

        if not mesh_motion_solver_settings.Has("mesh_velocity_calculation"):
            # add empty settings in case the user did not specify anything
            mesh_motion_solver_settings.AddEmptyValue("mesh_velocity_calculation")

        mesh_vel_calc_settings = mesh_motion_solver_settings["mesh_velocity_calculation"]

        fluid_solver_type = fluid_solver_settings["solver_type"].GetString()
        if fluid_solver_type == "monolithic" or fluid_solver_type == "Monolithic":
            if fluid_solver_settings.Has("time_scheme"):
                time_scheme_fluid = fluid_solver_settings["time_scheme"].GetString()
                alpha_fluid = fluid_solver_settings["alpha"].GetDouble()
                if mesh_vel_calc_settings.Has("time_scheme"):
                    time_scheme_mesh_vel = mesh_vel_calc_settings["time_scheme"].GetString()
                    if time_scheme_fluid != time_scheme_mesh_vel:
                        info_msg  = '"time_scheme" of the fluid (' + time_scheme_fluid
                        info_msg += ') is different from the\n"time_scheme" used for the '
                        info_msg += 'calculation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                        KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)
                else:
                    mesh_vel_calc_settings.AddValue("time_scheme", fluid_solver_settings["time_scheme"])
                    info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncalculation of the '
                    info_msg += 'mesh-velocity to "' + time_scheme_fluid + '" to be consistent with the\n'
                    info_msg += '"time_scheme" of the fluid'
                    KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)
                if not mesh_vel_calc_settings["time_scheme"].GetString().startswith("bdf"):
                    if mesh_vel_calc_settings.Has("alpha_m"):
                        alpha_mesh_vel = mesh_vel_calc_settings["alpha_m"].GetDouble()
                        if abs(alpha_fluid-alpha_mesh_vel) > 1e-12:
                            info_msg  = '"alpha" of the fluid (' + str(alpha_fluid)
                            info_msg += ') is different from the\n"alpha_m" used for the '
                            info_msg += 'calculation of the mesh-velocity (' + str(alpha_mesh_vel) + ')'
                            KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)
                    else:
                        mesh_vel_calc_settings.AddValue("alpha_m", fluid_solver_settings["alpha"])
                        info_msg  = 'setting "alpha_m" of the mesh-solver for the\ncalculation of the '
                        info_msg += 'mesh-velocity to "' + str(alpha_fluid) + '" to be consistent\nwith '
                        info_msg += '"alpha" of the fluid'
                        KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)

        elif fluid_solver_type == "fractional_step" or fluid_solver_type == "FractionalStep":
            # currently fractional step always uses BDF2
            if mesh_vel_calc_settings.Has("time_scheme"):
                time_scheme_mesh_vel = mesh_vel_calc_settings["time_scheme"].GetString()
                if time_scheme_mesh_vel != "bdf2":
                    info_msg  = '"time_scheme" of the fluid (bdf2) '
                    info_msg += 'is different from the\n"time_scheme" used for the '
                    info_msg += 'calculation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                    KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)
            else:
                mesh_vel_calc_settings.AddEmptyValue("time_scheme").SetString("bdf2")
                info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncalculation of the '
                info_msg += 'mesh-velocity to "bdf2" to be consistent with the\n'
                info_msg += '"time_scheme" of the fluid'
                KM.Logger.PrintInfo("::[NavierStokesAleFluidSolver]::", info_msg)

        if not mesh_vel_calc_settings.Has("time_scheme"):
            mesh_vel_calc_settings.AddEmptyValue("time_scheme").SetString("UNSPECIFIED")
            warn_msg  = 'unknown "solver_type" of the fluid-solver, therefore '
            warn_msg += 'no automatic selection\nof "time_scheme" for the calculation '
            warn_msg += 'of the mesh-velocity performed'
            KM.Logger.PrintWarning("::[NavierStokesAleFluidSolver]::", warn_msg)
