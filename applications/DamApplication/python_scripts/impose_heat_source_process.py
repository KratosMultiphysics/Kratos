import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

## This process sets the value of water loads.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeHeatSourceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeHeatSourceProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if ("NoorzaiHeatFlux2D" in settings["model_part_name"].GetString()) or ("NoorzaiHeatFlux3D" in settings["model_part_name"].GetString()):
            self.components_process_list.append(KratosDam.DamNoorzaiHeatFluxProcess(model_part, settings))

        if ("AzenhaHeatFlux2D" in settings["model_part_name"].GetString()) or ("AzenhaHeatFlux3D" in settings["model_part_name"].GetString()):
            self.components_process_list.append(KratosDam.DamAzenhaHeatFluxProcess(model_part, settings))

    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
