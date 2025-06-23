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

    class DifferentialOperator(KM.Process):
        '''Auxiliary class to execute the derivatives recovery.'''

        @staticmethod
        def GetDefaultParameters():
            return KM.Parameters("""
            {
                "operation"           : "gradient",
                "primitive_variable"  : "",
                "derivative_variable" : "",
                "buffer_step"         : 0,
                "process_step"        : "ExecuteFinalizeSolutionStep"
            }
            """)

        __operations = {
            "gradient"   : "RecoverGradient",
            "divergence" : "RecoverDivergence",
            "laplacian"  : "RecoverLaplacian"
        }

        def __init__(self, settings, recovery_tool, model_part):
            KM.Process.__init__(self)
            settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
            self.operation = self.__operations[settings["operation"].GetString()]
            self.primitive_variable = KM.KratosGlobals.GetVariable(settings["primitive_variable"].GetString())
            self.derivative_variable = KM.KratosGlobals.GetVariable(settings["derivative_variable"].GetString())
            self.buffer_step = settings["buffer_step"].GetInt()
            self.process_step = settings["process_step"].GetString()
            self.model_part = model_part
            self.differentiation = getattr(recovery_tool, self.operation)
            setattr(self, self.process_step, self.ExecuteDifferentiation)

        def ExecuteDifferentiation(self):
            self.differentiation(self.model_part, self.primitive_variable, self.derivative_variable, self.buffer_step)

        def Check(self):
            if self.operation == "RecoverGradient":
                if not isinstance(self.primitive_variable, KM.DoubleVariable):
                    raise Exception("The primitive variable of a gradient should be a scalar")
                if not isinstance(self.derivative_variable, KM.Array1DVariable3):
                    raise Exception("The derivative variable of a gradient should be a vector")
            if self.operation == "RecoverDivergence":
                if not isinstance(self.primitive_variable, KM.Array1DVariable3):
                    raise Exception("The primitive variable of a divergence should be a vector")
                if not isinstance(self.derivative_variable, KM.DoubleVariable):
                    raise Exception("The derivative variable of a divergence should be a scalar")
            if self.operation == "RecoverLaplacian":
                is_scalar = isinstance(self.primitive_variable, KM.DoubleVariable)
                is_vector = isinstance(self.primitive_variable, KM.Array1DVariable3)
                if is_scalar:
                    if not isinstance(self.derivative_variable, KM.DoubleVariable):
                        raise Exception("The derivative variable of a scalar laplacian should be a scalar")
                elif is_vector:
                    if not isinstance(self.primitive_variable, KM.Array1DVariable3):
                        raise Exception("The derivative variable of a vector laplacian should be a vector")
                else:
                    raise Exception("The primitive and derivative variables of a laplacian must be scalar or vector")


    def __init__(self, model, settings ):
        """Construct the DerivativesRecoveryProcess."""

        KM.Process.__init__(self)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.settings = settings
        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())

        domain_size = self.model_part.ProcessInfo.GetValue(KM.DOMAIN_SIZE)
        if domain_size == 2:
            self.recovery_tool = SW.DerivativesRecoveryUtility2D
        else:
            self.recovery_tool = SW.DerivativesRecoveryUtility3D

        self.operations = []
        for operation_settings in self.settings["list_of_operations"].values():
            operation = self.DifferentialOperator(operation_settings, self.recovery_tool, self.model_part)
            self.operations.append(operation)

    def Check(self):
        self.recovery_tool.Check(self.model_part)
        for operation in self.operations:
            operation.Check()

    def ExecuteInitialize(self):
        if self.settings["compute_neighbors"].GetBool():
            KM.FindGlobalNodalNeighboursProcess(self.model_part).Execute()
        self.recovery_tool.CalculatePolynomialWeights(self.model_part)

        for operation in self.operations:
            operation.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for operation in self.operations:
            operation.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        if self.settings["update_mesh_topology"].GetBool():
            self.ExecuteInitialize()

        for operation in self.operations:
            operation.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for operation in self.operations:
            operation.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        for operation in self.operations:
            operation.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        for operation in self.operations:
            operation.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        for operation in self.operations:
            operation.ExecuteFinalize()
