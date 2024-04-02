# Import Python modules
import json
import numpy
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return SaveROMCoefficientsProcess(model, settings["Parameters"])

class SaveROMCoefficientsProcess(KratosMultiphysics.OutputProcess):
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
            err_msg = "Unknown value \'{}\' for \'snapshots_control_type\'. Available options are \'time\' and \'step\'.".format(snapshots_control_type)
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


    @classmethod
    def GetDefaultParameters(self):
        default_settings = KratosMultiphysics.Parameters("""{
            "help": "A process to save the generalized coordinates (q) of the Reduced-Order Model (ROM) to disk.",
            "model_part_name": "",
            "snapshots_control_type": "step",
            "snapshots_interval": 1.0,
            "rom_coefficients_output_name": "RomCoefficientsSnapshots",
            "rom_coefficients_output_folder" : "rom_data"
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
        """Append the current ROM state's generalized coordinates (q) to the historical snapshots list.
        Note: 'current_rom_state' represents the generalized coordinates (q) at the current timestep."""
        current_rom_state = self.model_part.GetValue(KratosMultiphysics.RomApplication.ROM_SOLUTION_INCREMENT)
        self.rom_snapshots.append(current_rom_state)


        # Schedule next snapshot output
        if self.snapshots_interval > 0.0: # Note: if == 0, we'll just always print
            if self.snapshots_control_is_time:
                time = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
                while self.__GetPrettyFloat(self.next_output) <= time:
                    self.next_output += self.snapshots_interval
            else:
                step = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.STEP])
                while self.next_output <= step:
                    self.next_output += self.snapshots_interval


    def _PrintRomBasis(self):
        #Convert list to numpy
        self.rom_snapshots = numpy.array(self.rom_snapshots)

        # Create the folder if it doesn't already exist
        if not self.rom_coefficients_output_folder.exists():
            self.rom_coefficients_output_folder.mkdir(parents=True)

        numpy.save(self.rom_coefficients_output_folder / f"{self.rom_coefficients_output_name}.npy", self.rom_snapshots)

    def ExecuteFinalize(self):
        self._PrintRomBasis()

    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number