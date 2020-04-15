from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.MeshingApplication as MSH

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.eulerian_primitive_var_solver import EulerianPrimitiveVarSolver
from KratosMultiphysics.MeshingApplication.multiscale_refining_process import MultiscaleRefiningProcess

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
        self.condition_name = "LineCondition"
        self.min_buffer_size = 2

        # Initialize the multigrid process. It creates the model part
        self.multigrid = MultiscaleRefiningProcess(model, settings["multigrid_settings"])
        self.main_model_part = self.multigrid.GetRefinedModelPart()

        domain_size = self.settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)

        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def ImportModelPart(self):
        if self.main_model_part.ProcessInfo[MSH.SUBSCALE_INDEX] == 0:
            # Default implementation in the base class
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if self.main_model_part.ProcessInfo[MSH.SUBSCALE_INDEX] == 0:
            super(MultigridSolver, self).PrepareModelPart()
        self.multigrid.PrepareModelPart() # It creates the cpp utility instance

    def AdvanceInTime(self, current_time):
        divisions = 2**(self.GetComputingModelPart().GetValue(MSH.SUBSCALE_INDEX) * self.multigrid.number_of_divisions_at_subscale)
        dt = self._ComputeDeltaTime() / divisions
        new_time = current_time + dt

        self.GetComputingModelPart().CloneTimeStep(new_time)
        self.GetComputingModelPart().ProcessInfo[KM.STEP] += 1

        self.multigrid.ExecuteInitializeSolutionStep()

        KM.Logger.PrintInfo("::[Multigrid solver]::", "Subscale Index:", self.GetComputingModelPart().GetValue(MSH.SUBSCALE_INDEX))

        return new_time

    def GetComputingModelPart(self):
        return self.multigrid.GetRefinedModelPart()

    def InitializeSolutionStep(self):
        if self.GetComputingModelPart().NumberOfElements() != 0:
            super(MultigridSolver, self).InitializeSolutionStep()

    def Predict(self):
        if self.GetComputingModelPart().NumberOfElements() != 0:
            super(MultigridSolver, self).Predict()

    # def SolveSolutionStep(self):
    #     if self._TimeBufferIsInitialized():
    #         # if self.GetComputingModelPart().NumberOfElements() != 0:
    #         is_converged = self.solver.SolveSolutionStep()
    #         return is_converged

    def FinalizeSolutionStep(self):
        if self.GetComputingModelPart().NumberOfElements() != 0:
            super(MultigridSolver, self).FinalizeSolutionStep()

    def Finalize(self):
        self.multigrid.ExecuteFinalize()
