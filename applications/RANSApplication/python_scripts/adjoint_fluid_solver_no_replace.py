from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import applications
import KratosMultiphysics as Kratos

# Import base class file
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    from KratosMultiphysics.FluidDynamicsApplication.adjoint_vmsmonolithic_solver import AdjointVMSMonolithicSolver
else:
    msg = "RANSApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)


class AdjointFluidSolverNoReplace(AdjointVMSMonolithicSolver):
    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]:
            ## Set fluid properties from materials json file
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                Kratos.Logger.PrintWarning(
                    self.__class__.__name__,
                    "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json."
                )
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Model reading finished.")

    def GetComputingModelPart(self):
        return self.main_model_part
