import numpy as np
from pathlib import Path
import KratosMultiphysics
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object.")
    return ProjectedResidualsOutputProcess(model, settings["Parameters"])

class ProjectedResidualsOutputProcess(KratosMultiphysics.OutputProcess):
    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.settings = settings
        self.model_part = model[settings["model_part_name"].GetString()]
        self.output_path = Path(settings["output_path"].GetString())
        self.output_path.mkdir(parents=True, exist_ok=True)

        # Setup OutputController
        controller_settings = KratosMultiphysics.Parameters("""{}""")
        controller_settings.AddString("model_part_name", self.model_part.FullName())
        controller_settings.AddValue("output_control_type", settings["output_control_type"])
        controller_settings.AddValue("output_interval", settings["output_interval"])
        self.controller = KratosMultiphysics.OutputController(self.model_part.GetModel(), controller_settings)

        self.hrom_training_utility = None

    @classmethod
    def GetDefaultParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "model_part_name": "",
            "output_control_type": "step",
            "output_interval": 1,
            "output_path": "rom_data/Residuals",
            "sub_solver_name": "",
            "rom_basis_output_name": "RomParameters",
            "rom_basis_output_folder": "rom_data",
            "range_of_elements_to_export": [0, -1]
        }""")

    def SetDependencies(self, main_solver, rom_parameters):
        """Called by the Analysis Stage to inject the Python solver object."""
        sub_solver_name = self.settings["sub_solver_name"].GetString()

        if sub_solver_name != "":
            target_solver = getattr(main_solver, sub_solver_name)
        else:
            target_solver = main_solver

        # 1. Inject the parameters just like your analysis stage used to do
        if not rom_parameters.Has("rom_basis_output_name"):
            rom_parameters.AddString("rom_basis_output_name", self.settings["rom_basis_output_name"].GetString())
        else:
            rom_parameters["rom_basis_output_name"].SetString(self.settings["rom_basis_output_name"].GetString())

        if not rom_parameters.Has("rom_basis_output_folder"):
            rom_parameters.AddString("rom_basis_output_folder", self.settings["rom_basis_output_folder"].GetString())
        else:
            rom_parameters["rom_basis_output_folder"].SetString(self.settings["rom_basis_output_folder"].GetString())

        # 2. Instantiate the utility
        self.hrom_training_utility = HRomTrainingUtility(target_solver, rom_parameters)

    def IsOutputStep(self):
        return self.controller.Evaluate()

    def PrintOutput(self):
        if self.hrom_training_utility is not None:
            # 1. Extract the FULL math data
            res_mat = self.hrom_training_utility.GetCurrentResidualsProjected()

            # 2. Parse the export range directly from the process settings
            export_range = self.settings["range_of_elements_to_export"].GetVector()
            start = int(export_range[0])
            end = int(export_range[1]) if int(export_range[1]) >= 0 else None

            # Slice the array
            sliced_res_mat = res_mat[start:end, :]

            # 3. Determine file label
            if self.settings["output_control_type"].GetString() == "time":
                time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                file_label = f"{time:.7f}"
            else:
                step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
                file_label = f"{step}"

            # 4. Perform the I/O
            file_path = self.output_path / f"Residual_{file_label}.npy"
            np.save(file_path, sliced_res_mat)

        self.controller.Update()
