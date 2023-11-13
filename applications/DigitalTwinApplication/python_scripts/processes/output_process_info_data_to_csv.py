from pathlib import Path
import KratosMultiphysics as Kratos

def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if( not isinstance(settings, Kratos.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return OutputProcessInfoDataToCSV(model, settings["Parameters"])

class OutputProcessInfoDataToCSV(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        Kratos.OutputProcess.__init__(self)

        # Check default settings
        default_settings = Kratos.Parameters(r'''{
            "model_part_name"  : "",
            "list_of_variables": [],
            "output_prefix"    : ""
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.main_model_part = model[settings["model_part_name"].GetString()].GetRootModelPart()
        self.list_of_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in settings["list_of_variables"].GetStringArray()]
        self.output_prefix = Path(settings["output_prefix"].GetString())

    def PrintOutput(self):
        self.output_prefix.parent.mkdir(exist_ok=True, parents=True)

        process_info: Kratos.ProcessInfo = self.main_model_part.ProcessInfo
        step = process_info[Kratos.STEP]
        with open(str(self.output_prefix) + f"_{step:05d}.csv", "w") as file_output:
            # write the headers
            file_output.write(";".join([var.Name() for var in self.list_of_variables]))
            file_output.write("\n")
            file_output.write(";".join([str(process_info[var]) for var in self.list_of_variables]))
            file_output.write("\n")

