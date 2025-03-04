# Import Python modules
import json
import numpy
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM

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

        # Initialize variable for the model part
        self.model_part = None

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

        # Initialize first reduced snapshot
        self.init_rom_state = None

        # Initialize list of variables
        self.snapshot_variables_list = None

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
        "snapshot_solution_type": "incremental",  // Options: "total", "incremental"
        "nodal_unknowns": []
        }""")

        return default_settings
    

    def ConfigureOutputProcess(self, model_part, nodal_unknowns_list):

        # Get the model part from which the snapshots are to be retrieved
        self.model_part = model_part

        # Get list of variables that for the snapshot
        if len(nodal_unknowns_list) == 0:
            err_msg = "The snapshots matrix variables need to be specified by the user in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        if any(nodal_unknowns_list.count(var_name) > 1 for var_name in nodal_unknowns_list):
            err_msg = "There are repeated variables in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        nodal_unknowns_list.sort()

        self.snapshot_variables_list = []
        for var_name in nodal_unknowns_list:
            if not KratosMultiphysics.KratosGlobals.HasVariable(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not in KratosGlobals. Please check provided value.".format(var_name)
            if not KratosMultiphysics.KratosGlobals.GetVariableType(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(var_name)
            self.snapshot_variables_list.append(KratosMultiphysics.KratosGlobals.GetVariable(var_name))
    
    def ExecuteInitializeSolutionStep(self):
        if self.snapshot_solution_type=="incremental" and self.init_rom_state is None:                                      
            aux_init_snapshot = []
            for snapshot_var in self.snapshot_variables_list:
                aux_init_snapshot.append( numpy.array(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self.model_part.Nodes, snapshot_var, 0), copy=False ))
            init_snapshot = numpy.stack(aux_init_snapshot, axis=1).reshape(-1)
            self.init_rom_state = KratosROM.RomAuxiliaryUtilities.ProjectToReducedBasis(self.model_part, self.snapshot_variables_list, init_snapshot)

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

        aux_final_snapshot = []
        for snapshot_var in self.snapshot_variables_list:
            aux_final_snapshot.append( numpy.array(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self.model_part.Nodes, snapshot_var, 0), copy=False ))
        final_snapshot = numpy.stack(aux_final_snapshot, axis=1).reshape(-1)
        current_rom_state = KratosROM.RomAuxiliaryUtilities.ProjectToReducedBasis(self.model_part, self.snapshot_variables_list, final_snapshot)

        if self.snapshot_solution_type == "total":
            self.rom_snapshots.append(current_rom_state)
        elif self.snapshot_solution_type=="incremental":
            self.rom_snapshots.append(current_rom_state-self.init_rom_state)
            self.init_rom_state = current_rom_state

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