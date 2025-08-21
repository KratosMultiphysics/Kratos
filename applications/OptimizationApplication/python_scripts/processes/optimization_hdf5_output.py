import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
import KratosMultiphysics as Kratos

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationHDF5Output instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationHDF5Output(model, parameters["settings"], optimization_problem)

class OptimizationHDF5Output(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name": "MainModelPart",
            "file_settings": {
                "file_name": "<model_part_name>.h5",
                "file_access_mode": "read_write",
                "echo_level": 1
            },
            "element_data_value_settings": {
                "prefix": "/Step_<step>",
                "list_of_variables": []
            },
            "nodal_data_value_settings": {
                "prefix": "/Step_<step>",
                "list_of_variables" : []
            }
        }""")

        settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.parameters = settings
        self.optimization_problem = optimization_problem
        self.file_params = settings["file_settings"].Clone()

        self.model_part = model[settings["model_part_name"].GetString()]

    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        current_file_params = self.file_params.Clone()
        file_name = current_file_params["file_name"].GetString()
        file_name = file_name.replace("<model_part_name>", self.model_part.FullName())
        file_name = file_name.replace("<step>", str(self.optimization_problem.GetStep()))
        current_file_params["file_name"].SetString(file_name)
        with OpenHDF5File(current_file_params, self.model_part) as h5_file:
            current_settings = self.parameters["nodal_data_value_settings"].Clone()
            current_settings["prefix"].SetString(current_settings["prefix"].GetString().replace("<step>", str(self.optimization_problem.GetStep())))
            io = KratosHDF5.HDF5NodalDataValueIO(current_settings, h5_file)
            io.Write(self.model_part)

            current_settings = self.parameters["element_data_value_settings"].Clone()
            current_settings["prefix"].SetString(current_settings["prefix"].GetString().replace("<step>", str(self.optimization_problem.GetStep())))
            io = KratosHDF5.HDF5ElementDataValueIO(current_settings, h5_file)
            io.Write(self.model_part)


