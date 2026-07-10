import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

## This process sets the value of water loads.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeWaterLoadsConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeWaterLoadsConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if ("HydroLinePressure2D" in settings["model_part_name"].GetString()) or ("HydroSurfacePressure3D" in settings["model_part_name"].GetString()):

            self.components_process_list.append(KratosDam.DamHydroConditionLoadProcess(model_part, settings))

        if ("StraightUpliftLinePressure2D" in settings["model_part_name"].GetString()) or ("StraightUpliftSurfacePressure3D" in settings["model_part_name"].GetString()):

            joint_model_part = Model["MainModelPart.Parts_" + settings["joint_group_name"].GetString()]
            self.components_process_list.append(KratosDam.DamUpliftConditionLoadProcess(model_part, joint_model_part, settings))

        if "CircularUpliftSurfacePressure3D" in settings["model_part_name"].GetString():

            joint_model_part = Model["MainModelPart.Parts_" + settings["joint_group_name"].GetString()]
            self.components_process_list.append(KratosDam.DamUpliftCircularConditionLoadProcess(model_part, joint_model_part, settings))

        if ("HydroDynamicLinePressure2D" in settings["model_part_name"].GetString()) or ("HydroDynamicSurfacePressure3D" in settings["model_part_name"].GetString()):

            self.components_process_list.append(KratosDam.DamWestergaardConditionLoadProcess(model_part, settings))

    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()


