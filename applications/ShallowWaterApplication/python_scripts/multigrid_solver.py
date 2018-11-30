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
        self.multigrid.ExecuteInitialize()
        super(MultigridSolver, self).Initialize()

    # def ExecuteBeforeSolutionLoop(self):
    #     self.multigrid.ExecuteBeforeSolutionLoop()
    #     self.solver.Initialize()

    def ImportModelPart(self):
        print('------------------- THE CURRENT SUBSCALE INDEX -----------------------')
        print(self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX])
        if self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] == 0:
            # Default implementation in the base class
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] == 0:
            super(MultigridSolver, self).PrepareModelPart()
        self.multigrid.PrepareModelPart()

    def AdvanceInTime(self, current_time):
        divisions = self.main_model_part.ProcessInfo[Meshing.SUBSCALE_INDEX] * 2**self.multigrid.number_of_divisions_at_subscale
        dt = self._ComputeDeltaTime() / 2**divisions
        new_time = current_time + dt

        self._GetComputingModelPart().CloneTimeStep(new_time)
        self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1

        self.multigrid.ExecuteInitializeSolutionStep()

        return new_time

    def _GetComputingModelPart(self):
        return self.multigrid.GetRefinedModelPart()

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            if self._GetComputingModelPart().NumberOfElements() != 0:
                self.solver.InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            if self._GetComputingModelPart().NumberOfElements() != 0:
                self.solver.Predict()

    # def SolveSolutionStep(self):
    #     if self._TimeBufferIsInitialized():
    #         # if self._GetComputingModelPart().NumberOfElements() != 0:
    #         is_converged = self.solver.SolveSolutionStep()
    #         return is_converged

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            if self._GetComputingModelPart().NumberOfElements() != 0:
                self.solver.FinalizeSolutionStep()
