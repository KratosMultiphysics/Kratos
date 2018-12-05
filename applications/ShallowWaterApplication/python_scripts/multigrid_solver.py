from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow
import KratosMultiphysics.MeshingApplication as Meshing

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
from eulerian_primitive_var_solver import EulerianPrimitiveVarSolver
from multiscale_refining_process import MultiscaleRefiningProcess

def CreateSolver(model, custom_settings):
    return MultigridSolver(model, custom_settings)

class MultigridSolver(EulerianPrimitiveVarSolver):

    def __init__(self, model, settings):
        settings = self._ValidateSettings(settings)

        self.model = model      # TODO: inherit from PythonSolver and use super
        self.settings = settings
        self.echo_level = self.settings["echo_level"].GetInt()

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = "EulerPrimVarElement"
        self.condition_name = "Condition"
        self.min_buffer_size = 2

        # Initialize the multigrid process. It creates the model part
        self.multigrid = MultiscaleRefiningProcess(model, settings["multigrid_settings"])
        self.main_model_part = self.multigrid.GetRefinedModelPart()

        domain_size = self.settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def Initialize(self):
        super(MultigridSolver, self).Initialize()

        # TODO: Remove this temporary GiD output
        if self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] > 0:
            #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
            gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
            multifile = KratosMultiphysics.MultiFileFlag.SingleFile #MultipleFiles
            deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
            self.gid_io = KratosMultiphysics.GidIO("results_subscale",gid_mode,multifile,deformed_mesh_flag,write_conditions)

            mesh_name = 0.0
            self.gid_io.InitializeMesh( mesh_name );
            self.gid_io.WriteMesh(self.main_model_part.GetMesh());
            self.gid_io.FinalizeMesh()
            self.gid_io.InitializeResults(mesh_name,self.main_model_part.GetMesh())

    def ImportModelPart(self):
        if self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] == 0:
            # Default implementation in the base class
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] == 0:
            super(MultigridSolver, self).PrepareModelPart()
        self.multigrid.PrepareModelPart() # It creates the cpp utility instance

    def AdvanceInTime(self, current_time):
        divisions = 2**(self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] * self.multigrid.number_of_divisions_at_subscale)
        dt = self._ComputeDeltaTime() / divisions
        new_time = current_time + dt

        self._GetComputingModelPart().CloneTimeStep(new_time)
        self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1

        self.multigrid.ExecuteInitializeSolutionStep()

        return new_time

    def _GetComputingModelPart(self):
        return self.multigrid.GetRefinedModelPart()

    def InitializeSolutionStep(self):

        # for node in self.main_model_part.Nodes:
        #     node.Free(KratosMultiphysics.VELOCITY_X)
        #     node.Free(KratosMultiphysics.VELOCITY_Y)

        if self._GetComputingModelPart().NumberOfElements() != 0:
            super(MultigridSolver, self).InitializeSolutionStep()

        n_fixed_h = 0
        n_fixed_v1 = 0
        n_fixed_v2 = 0
        for node in self.main_model_part.Nodes:
            if node.IsFixed(Shallow.HEIGHT):
                n_fixed_h += 1
            if node.IsFixed(KratosMultiphysics.VELOCITY_X):
                n_fixed_v1 += 1
            if node.IsFixed(KratosMultiphysics.VELOCITY_Y):
                n_fixed_v2 += 1
                node.Set(KratosMultiphysics.BLOCKED)
        print("Number of nodes with fixed HEIGHT : ", n_fixed_h)
        print("Number of nodes with fixed VELOCITY_X : ", n_fixed_v1)
        print("Number of nodes with fixed VELOCITY_Y : ", n_fixed_v2)

        if self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] > 0:
            self.gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, self.main_model_part.Nodes, self.main_model_part.ProcessInfo[KratosMultiphysics.TIME], 0)
            self.gid_io.WriteNodalResults(Shallow.HEIGHT, self.main_model_part.Nodes, self.main_model_part.ProcessInfo[KratosMultiphysics.TIME], 0)
            self.gid_io.WriteNodalFlags(KratosMultiphysics.INTERFACE, "INTERFACE", self.main_model_part.Nodes, self.main_model_part.ProcessInfo[KratosMultiphysics.TIME])
            self.gid_io.WriteNodalFlags(KratosMultiphysics.BLOCKED, "BLOCKED", self.main_model_part.Nodes, self.main_model_part.ProcessInfo[KratosMultiphysics.TIME])
            self.gid_io.Flush()


    def Predict(self):
        if self._GetComputingModelPart().NumberOfElements() != 0:
            super(MultigridSolver, self).Predict()

    # def SolveSolutionStep(self):
    #     if self._TimeBufferIsInitialized():
    #         # if self._GetComputingModelPart().NumberOfElements() != 0:
    #         is_converged = self.solver.SolveSolutionStep()
    #         return is_converged

    def FinalizeSolutionStep(self):
        if self._GetComputingModelPart().NumberOfElements() != 0:
            super(MultigridSolver, self).FinalizeSolutionStep()
