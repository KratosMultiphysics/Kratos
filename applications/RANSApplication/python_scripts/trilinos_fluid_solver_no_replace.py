from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.trilinos_navier_stokes_solver_vmsmonolithic import TrilinosNavierStokesSolverMonolithic

class TrilinosFluidSolverNoReplace(TrilinosNavierStokesSolverMonolithic):
    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)
            self._SetPhysicalProperties()

        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)

        self.distributed_model_part_importer.CreateCommunicators()

        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.PrepareModelPart()

        KratosMultiphysics.Logger.PrintInfo("TrilinosFluidSolverNoReplace", "Model reading finished.")

    def GetComputingModelPart(self):
        return self.main_model_part


