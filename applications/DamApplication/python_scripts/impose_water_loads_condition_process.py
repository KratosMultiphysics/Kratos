from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of water loads.

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeWaterLoadsConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeWaterLoadsConditionProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if ("HydroLinePressure2D" in settings["model_part_name"].GetString()) or ("HydroSurfacePressure3D" in settings["model_part_name"].GetString()):

            self.components_process_list.append(DamHydroConditionLoadProcess(model_part, settings))

        if ("StraightUpliftLinePressure2D" in settings["model_part_name"].GetString()) or ("StraightUpliftSurfacePressure3D" in settings["model_part_name"].GetString()):

            self.components_process_list.append(DamUpliftConditionLoadProcess(model_part, settings))

        if "CircularUpliftSurfacePressure3D" in settings["model_part_name"].GetString():

            self.components_process_list.append(DamUpliftCircularConditionLoadProcess(model_part, settings))

        if ("HydroDynamicLinePressure2D" in settings["model_part_name"].GetString()) or ("HydroDynamicSurfacePressure3D" in settings["model_part_name"].GetString()):

            self.components_process_list.append(DamWestergaardConditionLoadProcess(model_part, settings))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()


