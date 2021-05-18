from . import standard_recoverer

class LagrangianMaterialAccelerationRecoverer(standard_recoverer.StandardMaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        standard_recoverer.StandardMaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)
    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.RecoverLagrangianAcceleration(self.model_part)
