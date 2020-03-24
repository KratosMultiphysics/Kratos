class ResponseFunctionInterface(object):
    """The base interface for response functions. Each response function
    is able to calculate its response value and gradient.
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
