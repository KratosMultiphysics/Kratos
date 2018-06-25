import KratosMultiphysics

class PythonProcess(object):
    def __init__(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
        
    def ExecuteInitialize(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def ExecuteBeforeSolutionLoop(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def ExecuteInitializeSolutionStep(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def ExecuteBeforeOutputStep(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def ExecuteAfterOutputStep(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def ExecuteFinalize(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")
    
    def Clear(self):
        KratosMultiphysics.Logger.PrintWarning("PythonProcess", "You are deriving your process from PythonProcess instead of Process. This is deprecated, please change")


def Factory(settings, Model):
    return PythonProcess()
