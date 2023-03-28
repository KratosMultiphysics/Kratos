import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

def Factory(settings, Model):
    if(not isinstance(settings,Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeInputTableNodalYoungModulusProcess(Model, settings["Parameters"])

class ImposeInputTableNodalYoungModulusProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
        initial_value = settings["initial_value"].GetDouble()
        settings.RemoveValue("min_value")
        settings.RemoveValue("max_value")

        # Checks if the parameter "input_file_name" exixts. If not, it create it and define it as "".
        # This may be changed in the future, when the Nodal Reference Temp (input file) process will be fixed.

        if settings.Has("input_file_name"):
            input_file_name = settings["input_file_name"].GetString()
        else:
            input_file_name = ""

        if ((input_file_name == "") or (input_file_name == "- No file") or (input_file_name == "- Add new file")):
            self.table = KratosMultiphysics.PiecewiseLinearTable()
        else:
            self.table = KratosMultiphysics.PiecewiseLinearTable()
            with open(input_file_name,'r') as file_name:
                for j, line in enumerate(file_name):
                    file_1 = line.split(" ")
                    if (len(file_1)) > 1:
                        self.table.AddRow(float(file_1[0]), float(file_1[1]))

        self.process = KratosDam.DamInputTableNodalYoungModulusProcess(model_part, self.table, settings)


    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
