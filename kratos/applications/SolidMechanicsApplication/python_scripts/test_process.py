from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
import python_process

#TODO: Josep Maria, this is just a test process. It shall be removed very soon, once the construction interface is agreed...

##all the processes python processes should be derived from "python_process"
class TestProcess(python_process.PythonProcess):
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
        print("inside solid_mechanics test process")
        
    def ExecuteFinalize(self):
        pass
    
    def Clear(self):
        pass


    
def Factory(settings, Model):
    return TestProcess()