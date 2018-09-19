from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")
import KratosMultiphysics.MeshMovingApplication

# Other imports
from python_solver import PythonSolver
import python_solvers_wrapper_mesh_motion


def CreateSolver(model, solver_settings, parallelism):
    return ALEFluidSolver(model, solver_settings, parallelism)


class ALEFluidSolver(PythonSolver):
    def __init__(self, model, solver_settings, parallelism):
        super(ALEFluidSolver, self).__init__(model, solver_settings)
        mesh_motion_solver_settings = solver_settings["ale_settings"].Clone()
        solver_settings.RemoveValue("ale_settings")

        fluid_model_part_name = solver_settings["model_part_name"].GetString()
        if not self.model.HasModelPart(fluid_model_part_name):
            fluid_mesh_model_part = KratosMultiphysics.ModelPart(fluid_model_part_name)
            self.model.AddModelPart(fluid_mesh_model_part)

        ## Creating the fluid solver
        self.fluid_solver = self._CreateFluidSolver(solver_settings, parallelism)

        ## Creating the mesh-motion solver
        if not mesh_motion_solver_settings.Has("echo_level"):
            mesh_motion_solver_settings.AddValue("echo_level", solver_settings["echo_level"])

        if mesh_motion_solver_settings.Has("model_part_name"):
            if not fluid_model_part_name == mesh_motion_solver_settings["model_part_name"].GetString():
                raise Exception('Fluid- and Mesh-Solver have to use the same "model_part_name"!')
        else:
            mesh_motion_solver_settings.AddValue("model_part_name", solver_settings["model_part_name"])

        domain_size = solver_settings["domain_size"].GetInt()
        if mesh_motion_solver_settings.Has("domain_size"):
            mesh_motion_domain_size = mesh_motion_solver_settings["domain_size"].GetInt()
            if not domain_size == mesh_motion_domain_size:
                raise Exception('Fluid- and Mesh-Solver have to use the same "domain_size"!')
        else:
            mesh_motion_solver_settings.AddValue("domain_size", solver_settings["domain_size"])

        self.ale_interface_parts_params = None
        if mesh_motion_solver_settings.Has("ale_interface_parts"):
            self.ale_interface_parts_params = mesh_motion_solver_settings["ale_interface_parts"].Clone()
            if not self.ale_interface_parts_params.IsArray():
                raise Exception('"ale_interface_parts" has to be provided as a list!')
            if self.ale_interface_parts_params.size() != domain_size:
                err_msg  = '"ale_interface_parts" has to have the same number of '
                err_msg += 'components as the domain_size (' + str(domain_size) + ')!'
                raise Exception(err_msg)
            mesh_motion_solver_settings.RemoveValue("ale_interface_parts")

        self.mesh_motion_solver = python_solvers_wrapper_mesh_motion.CreateSolverByParameters(
            model, mesh_motion_solver_settings, parallelism)

        # Getting the min_buffer_size from both solvers
        # and assigning it to the fluid_solver, bcs this one handles the model_part
        self.fluid_solver.min_buffer_size = max(self.fluid_solver.GetMinimumBufferSize(),
                                                self.mesh_motion_solver.GetMinimumBufferSize())

        self.is_printing_rank = self.fluid_solver._IsPrintingRank()

        # TODO move to "Check"?
        if (self.mesh_motion_solver.settings["calculate_mesh_velocities"].GetBool() == False
            and self.is_printing_rank):
            info_msg = "Mesh velocities are not being computed in the Mesh solver!"
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        if (self.fluid_solver.settings["compute_reactions"].GetBool() == False
            and self.is_printing_rank):
            info_msg = "Reactions are not being computed in the Fluid solver!"
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        # TODO once the different computations of the Mehs-Vel are implemented,
        # check if the time schemes are consistent (in fluid and for the computation
        # of the MESH_VELOCITY)

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Construction finished")

    def AddVariables(self):
        self.mesh_motion_solver.AddVariables()
        self.fluid_solver.AddVariables()
        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Variables Added")

    def AddDofs(self):
        self.mesh_motion_solver.AddDofs()
        self.fluid_solver.AddDofs()
        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "DOFs Added")

    def Initialize(self):
        # Saving the ALE-interface-parts for later and
        # fixing the MESH_DISPLACEMENT on the ALE-interface-parts
        # this has to be done AFTER reading the ModelPart
        domain_size = self.settings["domain_size"].GetInt()
        self.ale_interface_parts_by_components = [[] for x in range(domain_size)]
        if self.ale_interface_parts_params is not None:

            mesh_disp_components = [
                KratosMultiphysics.MESH_DISPLACEMENT_X,
                KratosMultiphysics.MESH_DISPLACEMENT_Y,
                KratosMultiphysics.MESH_DISPLACEMENT_Z
            ]

            main_model_part_name = self.settings["model_part_name"].GetString()

            for i_dir in range(domain_size):
                for i_name in range(self.ale_interface_parts_params[i_dir].size()):
                    # only submodelparts of the MainModelPart can be used!
                    sub_model_part_name = self.ale_interface_parts_params[i_dir][i_name].GetString()
                    full_model_part_name = main_model_part_name + "." + sub_model_part_name
                    ale_interface_model_part = self.model[full_model_part_name]
                    self.ale_interface_parts_by_components[i_dir].append(ale_interface_model_part)

                    # Fixing the components of the MESH_DISPLACEMENT
                    KratosMultiphysics.VariableUtils().ApplyFixity(
                        mesh_disp_components[i_dir], True, ale_interface_model_part.Nodes)

        self.mesh_motion_solver.Initialize()
        self.fluid_solver.Initialize()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Finished initialization")

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
        # This is necessary in order to preserve imposed velocities
        # this saves the values of the (historical) VELOCITY
        # as (non-historical) MESH_VELOCITY
        # This is done AFTER the processes have called ExecuteInitializeSolutionStep
        # (where ususally BCs are imposed)
        for ale_parts_dir in self.ale_interface_parts_by_components:
            for mp in ale_parts_dir:
                KratosMultiphysics.VariableUtils().SaveVectorVar(
                    KratosMultiphysics.VELOCITY,
                    KratosMultiphysics.MESH_VELOCITY,
                    mp.Nodes)

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

        # Copy the MESH_VELOCITY to the VELOCITY (ALE) on the interface
        mesh_vel_components = [
            KratosMultiphysics.MESH_VELOCITY_X,
            KratosMultiphysics.MESH_VELOCITY_Y,
            KratosMultiphysics.MESH_VELOCITY_Z
        ]
        vel_components = [
            KratosMultiphysics.VELOCITY_X,
            KratosMultiphysics.VELOCITY_Y,
            KratosMultiphysics.VELOCITY_Z
        ]

        for i_dir, ale_parts_dir in enumerate(self.ale_interface_parts_by_components):
            for mp in ale_parts_dir:
                # TODO implement this in C++
                for node in mp.Nodes:
                    # retrieve the previously stored velocity (see InitializeSolutionStep)
                    imposed_vel = node.GetValue(mesh_vel_components[i_dir])
                    mesh_vel = node.GetSolutionStepValue(mesh_vel_components[i_dir])
                    node.SetSolutionStepValue(vel_components[i_dir], imposed_vel + mesh_vel)

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
        It can be overridden to create different fluid solvers
        '''
        KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")
        import python_solvers_wrapper_fluid
        return python_solvers_wrapper_fluid.CreateSolverByParameters(
            self.model, solver_settings, parallelism)
