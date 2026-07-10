import numpy
from pathlib import Path
import KratosMultiphysics

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return NumpyOutputProcess(model, settings["Parameters"])


class NumpyOutputProcess(KratosMultiphysics.OutputProcess):
    """A process to save to disk numpy arrays with nodal solutions. No mesh information is stored"""

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

        # Get the variables list to be used to get the snapshots matrix information
        # Note that we sort the snapshot variables list alphabetically
        # This is required in order to establish a consensum for the possible visualization model part projections
        nodal_unknowns = settings["nodal_unknowns"].GetStringArray()
        if len(nodal_unknowns) == 0:
            err_msg = "The snapshots matrix variables need to be specified by the user in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        if any(nodal_unknowns.count(var_name) > 1 for var_name in nodal_unknowns):
            err_msg = "There are repeated variables in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        nodal_unknowns.sort()

        self.snapshot_variables_list = []
        for var_name in nodal_unknowns:
            if not KratosMultiphysics.KratosGlobals.HasVariable(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not in KratosGlobals. Please check provided value.".format(var_name)
            if not KratosMultiphysics.KratosGlobals.GetVariableType(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(var_name)
            self.snapshot_variables_list.append(KratosMultiphysics.KratosGlobals.GetVariable(var_name))

        # Setup the output folder
        self.output_path = Path(settings["output_path"].GetString())

        # Initialize output interval data
        self.next_output = 0.0


    @classmethod
    def GetDefaultParameters(self):
        default_settings = KratosMultiphysics.Parameters("""{
            "help": "A process to save the solution variables entered in 'nodal_unknowns'. Solutions' collection interlals are: 'step' or 'time'",
            "model_part_name": "",
            "snapshots_control_type": "step",
            "snapshots_interval": 1.0,
            "nodal_unknowns": [],
            "output_path": "numpy_output"
        }""")

        return default_settings

    # Flushing interval then writes to disk every n steps, or at the end of the simulation

    def IsOutputStep(self):
        if self.snapshots_control_is_time:
            time = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return time >= self.__GetPrettyFloat(self.next_output)
        else:
            step = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.STEP])
            return step >= self.next_output


    def PrintOutput(self):

        if not self.output_path.exists():
            self.output_path.mkdir(parents=True)

        step =  self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        time =  self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        aux_data_array = []
        for snapshot_var in self.snapshot_variables_list:
            aux_data_array.append( numpy.array(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self.model_part.Nodes, snapshot_var, 0), copy=False ))
        numpy.save(self.output_path / f"solution{step}.npy", numpy.stack(aux_data_array, axis=1).reshape(-1,1))

        # Schedule the next output
        current = time if self.snapshots_control_is_time else step
        while self.__GetPrettyFloat(self.next_output) <= self.__GetPrettyFloat(current):
            self.next_output += self.snapshots_interval


    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number