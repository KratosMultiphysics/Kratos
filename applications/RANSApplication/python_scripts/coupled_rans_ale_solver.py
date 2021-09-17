import KratosMultiphysics as Kratos
from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver

if CheckIfApplicationsAvailable("MeshMovingApplication"):
    from KratosMultiphysics.MeshMovingApplication.ale_fluid_solver import AleFluidSolver
else:
    raise Exception("In importing the CoupledRANSALESolver: The solver requires the MeshMovingApplication, but this application is not available.")

def CreateSolver(model, solver_settings, parallelism):
    return CoupledRANSALESolver(model, solver_settings, parallelism)


""" Currently only supports monolithic FSI """
class CoupledRANSALESolver(AleFluidSolver):
    def _CreateFluidSolver(self, solver_settings, parallelism):
        return CoupledRANSSolver(self.model, solver_settings)

    @classmethod
    def _ManipulateFluidSolverSettingsForReactionsComputation(cls, fluid_solver_settings):
        pass

    @classmethod
    def _ManipulateMeshMotionSolverSettingsForMeshVelocityComputation(cls, fluid_solver_settings, mesh_motion_solver_settings):
        # 1. Ensure that the MESH_VELOCITY is computed
        if mesh_motion_solver_settings.Has("calculate_mesh_velocity"):
            if not mesh_motion_solver_settings["calculate_mesh_velocity"].GetBool():
                mesh_motion_solver_settings["calculate_mesh_velocity"].SetBool(True)
                Kratos.Logger.PrintWarning(cls.__name__, 'Mesh velocity calculation was deactivated. Switching "calculate_mesh_velocity" on')
        else:
            mesh_motion_solver_settings.AddEmptyValue("calculate_mesh_velocity").SetBool(True)

        # 2. Selecting the time-integration for the MESH_VELOCITY to be consistent with the fluid time-integration.
        #    By now the parameters of the FluidSolver have been validated, which means
        #    that the time-integration method used by the fluid can be queried

        if not mesh_motion_solver_settings.Has("mesh_velocity_calculation"):
            # add empty settings in case the user did not specify anything
            mesh_motion_solver_settings.AddEmptyValue("mesh_velocity_calculation")

        mesh_vel_calc_settings = mesh_motion_solver_settings["mesh_velocity_calculation"]
        time_scheme_fluid = fluid_solver_settings["time_scheme_settings"]["scheme_type"].GetString()
        if mesh_vel_calc_settings.Has("time_scheme"):
            time_scheme_mesh_vel = mesh_vel_calc_settings["time_scheme"].GetString()
            if time_scheme_fluid != time_scheme_mesh_vel:
                info_msg  = '"time_scheme" of the fluid (' + time_scheme_fluid
                info_msg += ') is different from the\n"time_scheme" used for the '
                info_msg += 'calculation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                Kratos.Logger.PrintInfo(cls.__name__, info_msg)
        else:
            mesh_vel_calc_settings.AddValue("time_scheme", fluid_solver_settings["time_scheme_settings"]["scheme_type"])
            info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncalculation of the '
            info_msg += 'mesh-velocity to "' + time_scheme_fluid + '" to be consistent with the\n'
            info_msg += '"time_scheme" of the fluid'
            Kratos.Logger.PrintInfo(cls.__name__, info_msg)

        if time_scheme_fluid == "bossak":
            alpha_fluid = fluid_solver_settings["time_scheme_settings"]["alpha_bossak"].GetDouble()
            if mesh_vel_calc_settings.Has("alpha_m"):
                alpha_mesh_vel = mesh_vel_calc_settings["alpha_m"].GetDouble()
                if abs(alpha_fluid-alpha_mesh_vel) > 1e-12:
                    info_msg  = '"alpha" of the fluid (' + str(alpha_fluid)
                    info_msg += ') is different from the\n"alpha_m" used for the '
                    info_msg += 'calculation of the mesh-velocity (' + str(alpha_mesh_vel) + ')'
                    Kratos.Logger.PrintInfo(cls.__name__, info_msg)
            else:
                mesh_vel_calc_settings.AddValue("alpha_m", fluid_solver_settings["time_scheme_settings"]["alpha_bossak"])
                info_msg  = 'setting "alpha_m" of the mesh-solver for the\ncalculation of the '
                info_msg += 'mesh-velocity to "' + str(alpha_fluid) + '" to be consistent\nwith '
                info_msg += '"alpha" of the fluid'
                Kratos.Logger.PrintInfo(cls.__name__, info_msg)

