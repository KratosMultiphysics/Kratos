from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Other imports
from python_solver import PythonSolver
import python_solvers_wrapper_mesh_motion

class AleFluidSolver(PythonSolver):
    def __init__(self, model, solver_settings, parallelism):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"                 : "ale_fluid",
            "fluid_solver_settings"       : { },
            "mesh_motion_solver_settings" : { },
            "echo_level"                  : 0,
            "ale_boundary_parts"          : [ ]
        }""")

        # cannot recursively validate because validation of fluid- and
        # mesh-motion-settings is done in corresponding solvers
        solver_settings.ValidateAndAssignDefaults(default_settings)

        super(AleFluidSolver, self).__init__(model, solver_settings)

        fluid_solver_settings       = self.settings["fluid_solver_settings"]
        mesh_motion_solver_settings = self.settings["mesh_motion_solver_settings"]

        fluid_model_part_name = fluid_solver_settings["model_part_name"].GetString()
        if not self.model.HasModelPart(fluid_model_part_name):
            model.CreateModelPart(fluid_model_part_name)

        ## Checking if reactions are being computed in the fluid
        if fluid_solver_settings.Has("compute_reactions"):
            if fluid_solver_settings["compute_reactions"].GetBool() == False:
                fluid_solver_settings["compute_reactions"].SetBool(True)
                warn_msg  = '"compute_reactions" is switched off for the fluid-solver, '
                warn_msg += 'switching it on!'
                KratosMultiphysics.Logger.PrintWarning("::[AleFluidSolver]::", warn_msg)
        else:
            fluid_solver_settings.AddEmptyValue("compute_reactions").SetBool(True)
            info_msg = 'Setting "compute_reactions" to true for the fluid-solver'
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", info_msg)

        ## Creating the fluid solver
        self.fluid_solver = self._CreateFluidSolver(solver_settings, parallelism)
        self.is_printing_rank = self.fluid_solver._IsPrintingRank()

        ## Creating the mesh-motion solver
        if not mesh_motion_solver_settings.Has("echo_level"):
            mesh_motion_solver_settings.AddValue("echo_level", self.settings["echo_level"])

        # Making sure the settings are consistent btw fluid and mesh-motion
        if mesh_motion_solver_settings.Has("model_part_name"):
            if not fluid_model_part_name == mesh_motion_solver_settings["model_part_name"].GetString():
                raise Exception('Fluid- and Mesh-Solver have to use the same "model_part_name"!')
        else:
            mesh_motion_solver_settings.AddValue("model_part_name", fluid_solver_settings["model_part_name"])

        domain_size = fluid_solver_settings["domain_size"].GetInt()
        if mesh_motion_solver_settings.Has("domain_size"):
            mesh_motion_domain_size = mesh_motion_solver_settings["domain_size"].GetInt()
            if not domain_size == mesh_motion_domain_size:
                raise Exception('Fluid- and Mesh-Solver have to use the same "domain_size"!')
        else:
            mesh_motion_solver_settings.AddValue("domain_size", fluid_solver_settings["domain_size"])

        self.mesh_motion_solver = python_solvers_wrapper_mesh_motion.CreateSolverByParameters(
            model, mesh_motion_solver_settings, parallelism)

        # Getting the min_buffer_size from both solvers
        # and assigning it to the fluid_solver, bcs this one handles the model_part
        self.fluid_solver.min_buffer_size = max( [ self.fluid_solver.GetMinimumBufferSize(),
                                                   self.mesh_motion_solver.GetMinimumBufferSize(),
                                                   KratosMeshMoving.CalculateMeshVelocityUtility.
                                                   GetMinimumBufferSize(self.mesh_vel_comp_settings[
                                                       "time_scheme"].GetString()) ] )

        # TODO remove this once the MPI-logger is implemented
        self.is_printing_rank = self.fluid_solver._IsPrintingRank()

        # TODO once the different computations of the Mehs-Vel are implemented,
        # check if the time schemes are consistent (in fluid and for the computation
        # of the MESH_VELOCITY)
        # Then also the computation of the mesh-vel will be in a utility
        if (self.mesh_motion_solver.settings["calculate_mesh_velocities"].GetBool() == False
            and self.is_printing_rank):
            info_msg = "Mesh velocities are not being computed in the Mesh solver!"
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", info_msg)

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "Construction finished")

    def AddVariables(self):
        self.mesh_motion_solver.AddVariables()
        self.fluid_solver.AddVariables()

        # Adding Variables used for computation of Mesh-Velocity
        time_scheme = self.mesh_vel_comp_settings["time_scheme"].GetString()
        main_model_part = self.model[self.settings["model_part_name"].GetString()]
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        if not time_scheme.startswith("bdf"): # bdfx does not need MESH_ACCELERATION
            main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_ACCELERATION)

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "Variables Added")

    def AddDofs(self):
        self.mesh_motion_solver.AddDofs()
        self.fluid_solver.AddDofs()
        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "DOFs Added")

    def Initialize(self):
        # Saving the ALE-interface-parts for later
        # this can only be done AFTER reading the ModelPart
        self.ale_boundary_parts = []
        main_model_part_name = self.settings["fluid_solver_settings"]["model_part_name"].GetString()

        ale_boundary_parts_params = self.settings["ale_boundary_parts"]

        for i_name in range(ale_boundary_parts_params.size()):
            sub_model_part_name = ale_boundary_parts_params[i_name].GetString()
            full_model_part_name = main_model_part_name + "." + sub_model_part_name
            self.ale_boundary_parts.append(self.model[full_model_part_name])

        self.mesh_motion_solver.Initialize()
        self.fluid_solver.Initialize()

        self.calc_mesh_vel_util = KratosMeshMoving.CalculateMeshVelocityUtility(
            self.mesh_motion_solver.GetComputingModelPart(),
            self.mesh_vel_comp_settings)

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "Finished initialization")

    def ImportModelPart(self):
        self.fluid_solver.ImportModelPart() # only ONE solver imports the ModelPart

    def PrepareModelPart(self):
        # Doing it ONLY for the fluid solver (since this contains filling the buffer)
        self.fluid_solver.PrepareModelPart()

    def AdvanceInTime(self, current_time):
        # Doing it ONLY for the fluid solver
        return self.fluid_solver.AdvanceInTime(current_time)

    def Finalize(self):
        self.mesh_motion_solver.Finalize()
        self.fluid_solver.Finalize()

    def InitializeSolutionStep(self):
        self.mesh_motion_solver.InitializeSolutionStep()
        self.fluid_solver.InitializeSolutionStep()

    def Predict(self):
        self.mesh_motion_solver.Predict()
        self.fluid_solver.Predict()

    def FinalizeSolutionStep(self):
        self.mesh_motion_solver.FinalizeSolutionStep()
        self.fluid_solver.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        self.mesh_motion_solver.SolveSolutionStep()

        self.calc_mesh_vel_util.CalculateMeshVelocities()

        self._ApplyALEBoundaryCondition()

        self.fluid_solver.SolveSolutionStep()

    def Check(self):
        self.mesh_motion_solver.Check()
        self.fluid_solver.Check()

    def Clear(self):
        self.mesh_motion_solver.Clear()
        self.fluid_solver.Clear()

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart() # this is the same as the one used in the MeshSolver

    def GetFluidSolver(self):
        return self.fluid_solver

    def GetMeshMotionSolver(self):
        return self.mesh_motion_solver

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()


    def _CreateFluidSolver(self, solver_settings, parallelism):
        '''This function creates the fluid solver.
        It has to be overridden in derived classes
        '''
        raise Exception("Fluid solver creation must be implemented in the derived class.")


    def __ApplyALEBoundaryCondition(self):
        '''Copy the MESH_VELOCITY to the VELOCITY (ALE) on the ale-boundary
        '''
        for mp in self.ale_boundary_parts:
            KratosMultiphysics.VariableUtils().CopyVectorVar(
                KratosMultiphysics.MESH_VELOCITY,
                KratosMultiphysics.VELOCITY,
                mp.GetCommunicator().LocalMesh().Nodes)

    def _GetMeshVelocityComputationSettings(self, fluid_settings, mesh_motion_settings):
        # selecting the time-integration for the MESH_VELOCITY to be consistent
        # with the fluid time-integration
        # by now the parameters of the FluidSolver have been validated, which means
        # that the time-integration method used by the fluid can be queried

        mesh_vel_comp_settings = KratosMultiphysics.Parameters("""{ }""")

        # copying the parameters related to the mesh-vel computation
        if mesh_motion_settings.Has("time_scheme"):
            mesh_vel_comp_settings.AddValue("time_scheme", mesh_motion_settings["time_scheme"])
            mesh_motion_settings.RemoveValue("time_scheme")
        if mesh_motion_settings.Has("alpha_m"):
            mesh_vel_comp_settings.AddValue("alpha_m", mesh_motion_settings["alpha_m"])
            mesh_motion_settings.RemoveValue("alpha_m")
        if mesh_motion_settings.Has("alpha_f"):
            mesh_vel_comp_settings.AddValue("alpha_f", mesh_motion_settings["alpha_f"])
            mesh_motion_settings.RemoveValue("alpha_f")

        fluid_solver_type = fluid_settings["solver_type"].GetString()
        if fluid_solver_type == "monolithic" or fluid_solver_type == "Monolithic":
            if fluid_settings.Has("time_scheme"):
                time_scheme_fluid = fluid_settings["time_scheme"].GetString()
                alpha_fluid = fluid_settings["alpha"].GetDouble()
                if mesh_vel_comp_settings.Has("time_scheme"):
                    time_scheme_mesh_vel = mesh_vel_comp_settings["time_scheme"].GetString()
                    if time_scheme_fluid != time_scheme_mesh_vel and self.is_printing_rank:
                        info_msg  = '"time_scheme" of the fluid (' + time_scheme_fluid
                        info_msg += ') is different from the "time_scheme" used for the '
                        info_msg += 'computation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                else:
                    mesh_vel_comp_settings.AddValue("time_scheme", fluid_settings["time_scheme"])
                    if self.is_printing_rank:
                        info_msg  = 'setting "time_scheme" of the mesh-solver for the computation of the '
                        info_msg += 'mesh-velocity to "' + time_scheme_fluid + '" to be consistent with the '
                        info_msg += '"time_scheme" of the fluid'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                if mesh_vel_comp_settings.Has("alpha_m"):
                    alpha_mesh_vel = mesh_vel_comp_settings["alpha_"].GetDouble()
                    if abs(alpha_fluid-alpha_mesh_vel) > 1e-12 and self.is_printing_rank:
                        info_msg  = '"alpha" of the fluid (' + str(alpha_fluid)
                        info_msg += ') is different from the "alpha_m" used for the '
                        info_msg += 'computation of the mesh-velocity (' + str(alpha_mesh_vel) + ')'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                else:
                    mesh_vel_comp_settings.AddValue("alpha_m", fluid_settings["alpha"])
                    if self.is_printing_rank:
                        info_msg  = 'setting "alpha_m" of the mesh-solver for the computation of the '
                        info_msg += 'mesh-velocity to "' + str(alpha_fluid) + '" to be consistent with '
                        info_msg += '"alpha" of the fluid'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        elif fluid_solver_type == "fractional_step" or fluid_solver_type == "FractionalStep":
            # currently fractional step always uses BDF2
            if mesh_vel_comp_settings.Has("time_scheme"):
                time_scheme_mesh_vel = mesh_vel_comp_settings["time_scheme"].GetString()
                if time_scheme_mesh_vel != "bdf2" and self.is_printing_rank:
                    info_msg  = '"time_scheme" of the fluid (bdf2) '
                    info_msg += 'is different from the "time_scheme" used for the '
                    info_msg += 'computation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                    KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        if not mesh_vel_comp_settings.Has("time_scheme"):
            mesh_vel_comp_settings.AddEmptyValue("time_scheme").SetString("UNSPECIFIED")
            warn_msg  = 'unknown "solver_type" of the fluid-solver, therefore '
            warn_msg += 'no automatic selection\nof "time_scheme" for the computation '
            warn_msg += 'of the mesh-velocity performed'
            KratosMultiphysics.Logger.PrintWarning("::[ALEFluidSolver]::", warn_msg)

        return mesh_vel_comp_settings
