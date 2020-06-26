from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from . import standard_recoverer

class LagrangianMaterialAccelerationRecoverer(standard_recoverer.StandardMaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        standard_recoverer.StandardMaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)
    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.RecoverLagrangianAcceleration(self.model_part)
