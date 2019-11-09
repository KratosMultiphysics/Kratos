from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import applications
import KratosMultiphysics as Kratos

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic


class FluidSolverNoReplace(NavierStokesSolverMonolithic):
    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]:
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)
            self._SetPhysicalProperties()

        if not self.model.HasModelPart(
                self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)

        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.PrepareModelPart()

    def GetComputingModelPart(self):
        return self.main_model_part
