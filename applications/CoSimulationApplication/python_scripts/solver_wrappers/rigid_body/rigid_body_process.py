
import numpy as np


def Factory(solver, parameters):
    return _RigidBodyProcess(solver, parameters)


class _RigidBodyProcess():

    def __init__(self, solver, parameters):

        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        self.parameters = parameters
        self.solver = solver
    
    def Check(self):
        pass

    def ExecuteInitialize(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass
    
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    def ExecuteBeforeOutputStep(self):
        pass
    
    def ExecuteAfterOutputStep(self):
        pass
    
    def ExecuteFinalize(self):
        pass

    def GetDefaultParameters(self):
        pass
    
    def _AdaptTimeInterval(self, interval):
        if interval[1].IsString():
            if interval[1].GetString() == "End":
                interval[1].SetDouble(1e30)
        interval = np.array([interval[0].GetDouble(), interval[1].GetDouble()])
        return interval