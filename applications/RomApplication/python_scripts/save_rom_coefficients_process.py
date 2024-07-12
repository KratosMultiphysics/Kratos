# Import Python modules
import json
import numpy
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return SaveRomCoefficientsProcess(model, settings["Parameters"])

class SaveRomCoefficientsProcess(KratosMultiphysics.OutputProcess):
    """A process to save the generalized coordinates (q) of the Reduced-Order Model (ROM) to disk."""

    def __init__(self, model, settings):
        KratosMultiphysics.OutputProcess.__init__(self)

        # Validate input settings against defaults
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # Get the model part from which the snapshots are to be retrieved
        if not settings["model_part_name"].GetString():
            raise Exception("\'model_part_name\' not provided. Please specify the model part to get the snapshots from.")
        self.model_part = model[settings["model_part_name"].GetString()]

        # Set the snapshots output control and interval
        snapshots_control_type = settings["snapshots_control_type"].GetString()
        if snapshots_control_type == "time":
            self.snapshots_control_is_time = True
        elif snapshots_control_type == "step":
            self.snapshots_control_is_time = False
        else:
            err_msg = f"Unknown value \'{snapshots_control_type}\' for \'snapshots_control_type\'. Available options are \'time\' and \'step\'."
            raise Exception(err_msg)
        self.snapshots_interval = settings["snapshots_interval"].GetDouble()

        # Retrieve the folder name where the generalized coordinates will be written.
        self.rom_coefficients_output_folder = Path(settings["rom_coefficients_output_folder"].GetString())

        # Retrieve the numpy file name where the generalized coordinates will be written.
        self.rom_coefficients_output_name = Path(settings["rom_coefficients_output_name"].GetString())

        # Initialize output interval data
        self.next_output = 0.0

        # Initialize rom snapshots
        self.rom_snapshots = []

        # Initialize cumulative rom state
        self.cumulative_rom_state = None

        # Retrieve the user's preference for saving the ROM solution and ensure it's either 'total' or 'incremental'.
        self.snapshot_solution_type = settings["snapshot_solution_type"].GetString()
        if self.snapshot_solution_type not in ["total", "incremental"]:
            err_msg = "Unknown value '{}' for 'snapshot_solution_type'. Available options are 'total' and 'incremental'.".format(self.snapshot_solution_type)
            raise Exception(err_msg)



    @classmethod
    def GetDefaultParameters(self):
        default_settings = KratosMultiphysics.Parameters("""{
        "help": "This process saves the generalized coordinates (q) of the ROM. 'snapshot_solution_type' determines if the output is cumulative ('total') or stepwise ('incremental').",
        "model_part_name": "",
        "snapshots_control_type": "step",
        "snapshots_interval": 1.0,
        "rom_coefficients_output_folder": "rom_data",
        "rom_coefficients_output_name": "RomCoefficientsSnapshots",
        "snapshot_solution_type": "incremental"  // Options: "total", "incremental"
        }""")

        return default_settings

    def IsOutputStep(self):
        if self.snapshots_control_is_time:
            time = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return time >= self.__GetPrettyFloat(self.next_output)
        else:
            step = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.STEP])
            return step >= self.next_output

    def PrintOutput(self):
        """
        Append the current ROM state's generalized coordinates (q) to the historical snapshots list.
        If 'snapshot_solution_type' is set to 'total', the cumulative sum of q is saved, representing the total solution up to the current time step.
        If 'snapshot_solution_type' is 'incremental', only the solution increment (Δq) for the current time step is saved.
        """
        delta_q = self.model_part.GetValue(KratosMultiphysics.RomApplication.ROM_CURRENT_SOLUTION_TOTAL)
        if self.snapshot_solution_type == "total":
            if self.cumulative_rom_state is None:
                self.cumulative_rom_state = numpy.zeros_like(delta_q)
            self.cumulative_rom_state += delta_q
            self.rom_snapshots.append(numpy.copy(self.cumulative_rom_state))
        elif self.snapshot_solution_type=="incremental":
            self.rom_snapshots.append(delta_q)

        # Schedule the next output
        current = self.model_part.ProcessInfo[KratosMultiphysics.TIME if self.snapshots_control_is_time else KratosMultiphysics.STEP]
        while self.__GetPrettyFloat(self.next_output) <= self.__GetPrettyFloat(current):
            self.next_output += self.snapshots_interval

    def _PrintRomCoefficients(self):
        #Convert list to numpy
        self.rom_snapshots = numpy.array(self.rom_snapshots)

        # Create the folder if it doesn't already exist
        if not self.rom_coefficients_output_folder.exists():
            self.rom_coefficients_output_folder.mkdir(parents=True)

        numpy.save(self.rom_coefficients_output_folder / f"{self.rom_coefficients_output_name}.npy", self.rom_snapshots)

    def ExecuteFinalize(self):
        self._PrintRomCoefficients()

    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number