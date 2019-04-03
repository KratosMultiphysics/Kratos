from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of water loads.

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeHeatSourceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeHeatSourceProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if ("NoorzaiHeatFlux2D" in settings["model_part_name"].GetString()) or ("NoorzaiHeatFlux3D" in settings["model_part_name"].GetString()):
            self.components_process_list.append(DamNoorzaiHeatFluxProcess(model_part, settings))

        if ("AzenhaHeatFlux2D" in settings["model_part_name"].GetString()) or ("AzenhaHeatFlux3D" in settings["model_part_name"].GetString()):
            self.components_process_list.append(DamAzenhaHeatFluxProcess(model_part, settings))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
