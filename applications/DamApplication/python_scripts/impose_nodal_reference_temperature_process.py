from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *



def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeNodalReferenceTemperatureProcess(Model, settings["Parameters"])

class ImposeNodalReferenceTemperatureProcess(Process):
    
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
        initial_value = settings["initial_value"].GetDouble()
        
        if(initial_value == 0):
            input_file_name = settings["input_file_name"].GetString()
            self.table = PiecewiseLinearTable()
            with open(input_file_name,'r') as file_name:
                for j, line in enumerate(file_name):
                    file_1 = line.split(" ")
                    if (len(file_1)) > 1: 
                        self.table.AddRow(float(file_1[0]), float(file_1[1]))
        else:
            self.table = PiecewiseLinearTable()

        self.process = DamNodalReferenceTemperatureProcess(model_part, self.table, settings) 

                 
    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
