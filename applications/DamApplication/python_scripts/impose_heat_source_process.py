from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of water loads.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeHeatSourceProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeHeatSourceProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]
        
        self.components_process_list = []
        
        if "Noorzai" in settings["model_part_name"].GetString():
            self.components_process_list.append(DamNoorzaiHeatFluxProcess(model_part, settings))
                       
        if "Azenha" in settings["model_part_name"].GetString(): 
            self.components_process_list.append(DamAzenhaHeatFluxProcess(model_part, settings))
          
    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()