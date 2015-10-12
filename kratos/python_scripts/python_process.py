import KratosMultiphysics

class PythonProcess(object):
    def __init__(self):
        pass
        
    def ExecuteInitialize(self):
        pass
    
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass
    
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    def ExecuteBeforeOutputStep(self):
        pass
    
    def ExecuteAfterOutputStep(self):
        print("aaa")
        pass
    
    def ExecuteFinalize(self):
        pass
    
    def Clear(self):
        pass


def Factory(settings, Model):
    return PythonProcess()