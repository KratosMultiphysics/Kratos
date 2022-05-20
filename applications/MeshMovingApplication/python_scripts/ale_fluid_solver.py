# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from KratosMultiphysics.python_solver import PythonSolver, PhysicalSolver
import KratosMultiphysics.MeshMovingApplication.python_solvers_wrapper_mesh_motion as mesh_mothion_solvers_wrapper

class AleFluidSolver(PythonSolver):
    
    def __init__(self, model, solver_settings, parallelism):

        self._validate_settings_in_baseclass=True # To be removed eventually
        
        super().__init__(model, solver_settings)

        self.start_fluid_solution_time = self.settings["start_fluid_solution_time"].GetDouble()

        print("Assing Parallelism type for", self)
        self.parallelism = parallelism

        fluid_solver_settings       = self.settings["fluid_solver_settings"]
        mesh_motion_solver_settings = self.settings["mesh_motion_solver_settings"]

        fluid_model_part_name = fluid_solver_settings["model_part_name"].GetString()
        if not self.model.HasModelPart(fluid_model_part_name):
            model.CreateModelPart(fluid_model_part_name)

        # Derived class decides if the reactions should be computed or not
        self._ManipulateFluidSolverSettingsForReactionsComputation(self.settings["fluid_solver_settings"])

        # # Creating the fluid solver (DOWN)
        # self.GetPhysicalSolver("Fluid").Create(self.GetPhysicalSolver("Fluid").settings, parallelism)
        
        # Creating the mesh-motion solver
        if not self.settings["mesh_motion_solver_settings"].Has("echo_level"):
            self.settings["mesh_motion_solver_settings"].AddValue("echo_level", self.settings["echo_level"])

        # Making sure the settings are consistent btw fluid and mesh-motion
        if self.settings["mesh_motion_solver_settings"].Has("model_part_name"):
            if not fluid_model_part_name == self.settings["mesh_motion_solver_settings"]["model_part_name"].GetString():
                err_msg =  'Fluid- and Mesh-Solver have to use the same MainModelPart ("model_part_name")!\n'
                err_msg += 'Use "mesh_motion_parts" for specifying mesh-motion on sub-model-parts'
                raise Exception(err_msg)
        else:
            self.settings["mesh_motion_solver_settings"].AddValue("model_part_name", self.settings["fluid_solver_settings"]["model_part_name"])

        domain_size = self.settings["fluid_solver_settings"]["domain_size"].GetInt()
        if self.settings["mesh_motion_solver_settings"].Has("domain_size"):
            mesh_motion_domain_size = self.settings["mesh_motion_solver_settings"]["domain_size"].GetInt()
            if not domain_size == mesh_motion_domain_size:
                raise Exception('Fluid- and Mesh-Solver have to use the same "domain_size"!')
        else:
            self.settings["mesh_motion_solver_settings"].AddValue("domain_size", self.settings["fluid_solver_settings"]["domain_size"])

        # Derived class decides if the mesh velocities should be computed or not
        self._ManipulateMeshMotionSolverSettingsForMeshVelocityComputation(self.settings["fluid_solver_settings"], self.settings["mesh_motion_solver_settings"])

        # # Constructing the mesh-solver with the entire mesh (DOWN)
        # # if no submodelparts are specified then this is used for the computation of the mesh-motion
        # # otherwise it only adds the dofs and the variables (to the entire ModelPart!)
        # self.mesh_motion_solver_full_mesh = self.GetPhysicalSolver("MeshMotion").Create(
        #     model, self.GetPhysicalSolver("MeshMotion").settings, parallelism
        # )

        # # Getting the min_buffer_size from both solvers (DOWN)
        # # and assigning it to the fluid_solver, bcs this one handles the model_part
        # self.GetPhysicalSolver("Fluid")().min_buffer_size = max(
        #     self.GetPhysicalSolver("Fluid")().GetMinimumBufferSize(),
        #     self.mesh_motion_solver_full_mesh.GetMinimumBufferSize())

        print("CREATING SOLVERS")
        # Create the solvers with the arguments and the functions registered
        # for solver in self.GetPhysicalSolvers():
        #     solver.CreateSolver(*solver.arguments)

        KM.Logger.PrintInfo("::[AleFluidSolver]::", "Construction finished")

    def _RegisterPhysicalSolvers(self):
        super()._RegisterPhysicalSolvers()
        self.AddPhysicalSolver("MeshMotion", PhysicalSolver(
            self.settings["mesh_motion_solver_settings"],
            mesh_mothion_solvers_wrapper.CreateSolverByParameters,
            [self.model, self.settings["mesh_motion_solver_settings"], self.parallelism]
        ))


    def _CreatePhysicalSolvers(self):
        super()._CreatePhysicalSolvers()
        self.GetPhysicalSolver("Fluid")().min_buffer_size = max(
            self.GetPhysicalSolver("Fluid")().GetMinimumBufferSize(),
            self.GetPhysicalSolver("MeshMotion")().GetMinimumBufferSize()
        )

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver_type"                 : "ale_fluid",
            "start_fluid_solution_time"   : 0.0,
            "ale_boundary_parts"          : [ ],
            "mesh_motion_parts"           : [ ],
            "fluid_solver_settings"       : { },
            "mesh_motion_solver_settings" : { }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        self.GetPhysicalSolver("MeshMotion")().AddVariables()
        self.GetPhysicalSolver("Fluid")().AddVariables()

        KM.Logger.PrintInfo("::[AleFluidSolver]::", "Variables Added")

    def AddDofs(self):
        self.GetPhysicalSolver("MeshMotion")().AddDofs()
        self.GetPhysicalSolver("Fluid")().AddDofs()
        KM.Logger.PrintInfo("::[AleFluidSolver]::", "DOFs Added")

    def Initialize(self):
        # Getting the min_buffer_size from both solvers
        # and assigning it to the fluid_solver, bcs this one handles the model_part
        self.GetPhysicalSolver("Fluid")().min_buffer_size = max(
            self.GetPhysicalSolver("Fluid")().GetMinimumBufferSize(),
            self.GetPhysicalSolver("MeshMotion")().GetMinimumBufferSize())

        # Saving the ALE-interface-parts for later
        # this can only be done AFTER reading the ModelPart
        main_model_part_name = self.settings["fluid_solver_settings"]["model_part_name"].GetString()

        ale_boundary_parts_params = self.settings["ale_boundary_parts"]
        self.ale_boundary_parts = []
        for i_name in range(ale_boundary_parts_params.size()):
            sub_model_part_name = ale_boundary_parts_params[i_name].GetString()
            full_model_part_name = main_model_part_name + "." + sub_model_part_name
            self.ale_boundary_parts.append(self.model[full_model_part_name])

        mesh_motion_parts_params = self.settings["mesh_motion_parts"]
        self.mesh_motion_solvers = []
        if mesh_motion_parts_params.size() == 0:
            # the entire Fluid-ModelPart is used in the Mesh-Solver
            self.mesh_motion_solvers.append(self.GetPhysicalSolver("MeshMotion")())
        else:
            # SubModelParts of the Fluid-ModelPart are used in the Mesh-Solver
            # each SubModelPart has its own mesh-solver
            # Note that these solvers do NOT need to call AddVariables and AddDofs
            # since this is done already for the MainModelPart
            for i_name in range(mesh_motion_parts_params.size()):
                sub_model_part_name = mesh_motion_parts_params[i_name].GetString()
                if sub_model_part_name == main_model_part_name:
                    err_msg =  'The MainModelPart cannot be used as one of the Sub-Mesh-Solvers!\n'
                    err_msg += 'Remove "mesh_motion_parts" for specifying mesh-motion on the MainModelPart'
                    raise Exception(err_msg)
                full_model_part_name = main_model_part_name + "." + sub_model_part_name
                sub_mesh_solver_settings = self.settings["mesh_motion_solver_settings"].Clone()
                sub_mesh_solver_settings["model_part_name"].SetString(full_model_part_name)

                self.mesh_motion_solvers.append(mesh_mothion_solvers_wrapper.CreateSolverByParameters(
                    self.model, sub_mesh_solver_settings, self.parallelism))

        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Initialize()
        self.GetPhysicalSolver("Fluid")().Initialize()

        KM.Logger.PrintInfo("::[AleFluidSolver]::", "Finished initialization")

    def ImportModelPart(self):
        self.GetPhysicalSolver("Fluid")().ImportModelPart() # only ONE mesh_solver imports the ModelPart

    def PrepareModelPart(self):
        # Doing it ONLY for the fluid solver (since this contains filling the buffer)
        self.GetPhysicalSolver("Fluid")().PrepareModelPart()

    def AdvanceInTime(self, current_time):
        # Doing it ONLY for the fluid solver
        return self.GetPhysicalSolver("Fluid")().AdvanceInTime(current_time)

    def Finalize(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Finalize()
        self.GetPhysicalSolver("Fluid")().Finalize()

    def InitializeSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.InitializeSolutionStep()
        self.GetPhysicalSolver("Fluid")().InitializeSolutionStep()

    def Predict(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Predict()
        self.GetPhysicalSolver("Fluid")().Predict()

    def FinalizeSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.FinalizeSolutionStep()
        self.GetPhysicalSolver("Fluid")().FinalizeSolutionStep()

    def SolveSolutionStep(self):
        is_converged = True
        for mesh_solver in self.mesh_motion_solvers:
            is_converged &= mesh_solver.SolveSolutionStep()

        if self.GetPhysicalSolver("Fluid")().GetComputingModelPart().ProcessInfo[KM.TIME] >= self.start_fluid_solution_time:
            self.__ApplyALEBoundaryCondition()
            is_converged &= self.GetPhysicalSolver("Fluid")().SolveSolutionStep()

        return is_converged

    def Check(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Check()
        self.GetPhysicalSolver("Fluid")().Check()

    def Clear(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Clear()
        self.GetPhysicalSolver("Fluid")().Clear()

    def GetComputingModelPart(self):
        return self.GetPhysicalSolver("Fluid")().GetComputingModelPart() # this is the same as the one used in the MeshSolver

    def GetFluidSolver(self):
        return self.GetPhysicalSolver("Fluid")()

    def GetMeshMotionSolver(self):
        if len(self.mesh_motion_solvers) > 1:
            raise Exception('More than one mesh-motion-solver exists, please use "GetMeshMotionSolvers"')
        return self.mesh_motion_solvers[0]

    def GetMeshMotionSolvers(self):
        return self.mesh_motion_solvers

    def MoveMesh(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.MoveMesh()

    def _CreateFluidSolver(self, solver_settings, parallelism):
        '''This function creates the fluid solver.
        It has to be overridden in derived classes
        '''
        raise NotImplementedError("Fluid solver creation must be implemented in the derived class.")

    def __ApplyALEBoundaryCondition(self):
        '''Copy the MESH_VELOCITY to the VELOCITY (ALE) on the ale-boundary
        '''
        for mp in self.ale_boundary_parts:
            KM.VariableUtils().CopyVectorVar(
                KM.MESH_VELOCITY,
                KM.VELOCITY,
                mp.GetCommunicator().LocalMesh().Nodes)
            mp.GetCommunicator().SynchronizeVariable(KM.VELOCITY)

    @classmethod
    def _ManipulateFluidSolverSettingsForReactionsComputation(cls, fluid_solver_settings):
        raise NotImplementedError('"_ManipulateFluidSolverSettingsForReactionsComputation" has to be implemented in the derived class!')

    @classmethod
    def _ManipulateMeshMotionSolverSettingsForMeshVelocityComputation(cls, fluid_solver_settings, mesh_motion_solver_settings):
        raise NotImplementedError('"_ManipulateMeshMotionSolverSettingsForMeshVelocityComputation" has to be implemented in the derived class!')
