import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DerivativesRecoveryProcess(Model, settings["Parameters"])


class DerivativesRecoveryProcess(KM.Process):

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""
        {
            "model_part_name"          : "model_part",
            "list_of_operations"       : [],
            "compute_neighbors"        : true,
            "update_mesh_topology"     : false
        }
        """)

    class DifferentialOperator():
        '''Auxiliary class to execute the derivatives recovery.'''

        @staticmethod
        def GetDefaultParameters():
            return KM.Parameters("""
            {
                "operation"           : "gradient",
                "primitive_variable"  : "",
                "derivative_variable" : "",
                "buffer_step"         : 0
            }
            """)

        __operations = {
            "gradient"   : "RecoverGradient",
            "divergence" : "RecoverDivergence",
            "laplacian"  : "RecoverLaplacian"
        }

        def __init__(self, settings):
            settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
            self.operation = self.__operations[settings["operation"].GetString()]
            self.primitive_variable = KM.KratosGlobals.GetVariable(settings["primitive_variable"].GetString())
            self.derivative_variable = KM.KratosGlobals.GetVariable(settings["derivative_variable"].GetString())
            self.buffer_step = settings["buffer_step"].GetInt()

        def __call__(self, recovery_tool):
            getattr(recovery_tool, self.operation)(self.primitive_variable, self.derivative_variable, self.buffer_step)

        def Check(self):
            if self.operation is "RecoverGradient":
                if not isinstance(self.primitive_variable, KM.DoubleVariable):
                    raise Exception("The primitive variable of a gradient should be scalar type")
                if not isinstance(self.derivative_variable, KM.Array1DVariable3):
                    raise Exception("The derivative variable of a gradient should be vector type")



    def __init__(self, model, settings ):
        """Construct the DerivativesRecoveryProcess."""

        KM.Process.__init__(self)

        settings.ValidateAndAssignDefaults()
        self.settings = settings
        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())

        domain_size = self.model_part.ProcessInfo.GetValue(KM.DOMAIN_SIZE)
        if domain_size == 2:
            self.utility = SW.DerivativesRecoveryUtility2D
        else:
            self.utility = SW.DerivativesRecoveryUtility3D

        for operation_settings in self.settings["list_of_operations"]:
            operation_settings.ValidateAndAssignDefaults(self.GetDefaultOperationsParameters())