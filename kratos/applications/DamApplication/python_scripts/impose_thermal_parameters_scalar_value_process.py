from KratosMultiphysics import *

## This proces sets the value of a scalar variable using the ApplyConstantScalarValueProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeThemalParametersScalarValueProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeThemalParametersScalarValueProcess(ApplyConstantVectorValueProcess):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]
        
        thermal_density = settings["ThermalDensity"].GetDouble()
        conductivity = settings["Conductivity"].GetDouble()
        specific_heat = settings["SpecificHeat"].GetDouble()
        
        if thermal_density > 0.00001:
            # TODO: We have to assign the variable inside
            ApplyConstantScalarValueProcess.__init__(self,model_part, settings)
            print (thermal_density)
        if conductivity > 0.00001: 
            print (conductivity)
        if specific_heat > 0.00001:
            print (specific_heat)  
