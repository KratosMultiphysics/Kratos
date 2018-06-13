from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *



def Factory(settings, Model):
    if(not isinstance(settings,Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeNodalReferenceTemperatureProcess(Model, settings["Parameters"])

class ImposeNodalReferenceTemperatureProcess(Process):
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
        initial_value = settings["initial_value"].GetDouble()

        # Checks if the parameter "input_file_name" exixts. If not, it create it and define it as "".
        # This may be changed in the future, when the Nodal Reference Temp (input file) process will be fixed.

        if settings.Has("input_file_name"):
            input_file_name = settings["input_file_name"].GetString()
        else:
            input_file_name = ""

        if ((input_file_name == "") or (input_file_name == "- No file") or (input_file_name == "- Add new file")):
            self.table = PiecewiseLinearTable()
        else:
            self.table = PiecewiseLinearTable()
            with open(input_file_name,'r') as file_name:
                for j, line in enumerate(file_name):
                    file_1 = line.split(" ")
                    if (len(file_1)) > 1:
                        self.table.AddRow(float(file_1[0]), float(file_1[1]))

        self.process = DamNodalReferenceTemperatureProcess(model_part, self.table, settings)


    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
