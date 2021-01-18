from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This is process is related to add_mass processes

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SpecialConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class SpecialConditionProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
