class ResponseFunctionInterface(object):
    """
    This response function interface provides a unique interface for all possible ways
    to calculate the value and gradient of a response.

    The interface is designed to be used in e.g. optimization, where the value and gradient
    of a response is required, however the exact method of gradient calculation is of
    secondary importance.

    This might be done using e.g. adjoint sensitivity analysis capabilities of Kratos,
    or even a simple finite differencing method.

    (Do not confuse this class with the kratos/response_functions/adjoint_response_function.h,
    which is an implementation detail for the adjoint sensitivity analysis in Kratos)
    """

    def RunCalculation(self, calculate_gradient):
        self.Initialize()
        self.InitializeSolutionStep()
        self.CalculateValue()
        if calculate_gradient:
            self.CalculateGradient()
        self.FinalizeSolutionStep()
        self.Finalize()

    def Initialize(self):
        pass

    def UpdateDesign(self, updated_model_part, variable):
        pass

    def InitializeSolutionStep(self):
        pass

    def CalculateValue(self):
        raise NotImplementedError("CalculateValue needs to be implemented by the derived class")

    def CalculateGradient(self):
        raise NotImplementedError("CalculateGradient needs to be implemented by the derived class")

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    def GetValue(self):
        raise NotImplementedError("GetValue needs to be implemented by the derived class")

    def GetNodalGradient(self, variable):
        raise NotImplementedError("GetNodalGradient needs to be implemented by the derived class")

    def GetElementalGradient(self, variable):
        raise NotImplementedError("GetElementalGradient needs to be implemented by the derived class")

    def GetConditionalGradient(self, variable):
        raise NotImplementedError("GetConditionalGradient needs to be implemented by the derived class")

    def IsEvaluatedInFolder(self):
        return False
