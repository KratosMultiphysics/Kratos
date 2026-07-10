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

        # Get the model part from which the solutions are to be retrieved
        if not settings["model_part_name"].GetString():
            raise Exception("\'model_part_name\' not provided. Please specify the model part to get the solution from.")
        self.model_part = model[settings["model_part_name"].GetString()]

        # Set the solution output control and interval
        output_control_type = settings["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            err_msg = "Unknown value \'{}\' for \'output_control_type\'. Available options are \'time\' and \'step\'.".format(output_control_type)
            raise Exception(err_msg)
        self.output_interval = settings["output_interval"].GetDouble()

        nodal_results = settings["nodal_results"].GetStringArray()
        if len(nodal_results) == 0:
            err_msg = "The snapshots matrix variables need to be specified by the user in the \'nodal_results\' string array."
            raise Exception(err_msg)
        if any(nodal_results.count(var_name) > 1 for var_name in nodal_results):
            err_msg = "There are repeated variables in the \'nodal_results\' string array."
            raise Exception(err_msg)
        # The snapshot variables list is sorted alphabetically
        nodal_results.sort()

        self.snapshot_variables_list = []
        for var_name in nodal_results:
            if not KratosMultiphysics.KratosGlobals.HasVariable(var_name):
                err_msg = "\'{}\' variable in \'nodal_results\' is not in KratosGlobals. Please check provided value.".format(var_name)
                raise Exception(err_msg)
            if KratosMultiphysics.KratosGlobals.GetVariableType(var_name) != 'Double':
                err_msg = "\'{}\' variable in \'nodal_results\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(var_name)
                raise Exception(err_msg)
            self.snapshot_variables_list.append(KratosMultiphysics.KratosGlobals.GetVariable(var_name))

        # Setup the output folder
        self.output_path = Path(settings["output_path"].GetString())

        # Initialize output interval data
        self.next_output = 0.0


    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""{
            "help": "A process to save the solution variables entered in 'nodal_results'. output_control_type can be: 'step' or 'time' ",
            "model_part_name": "",
            "output_control_type": "step",
            "output_interval": 1.0,
            "nodal_results": [],
            "output_path": "numpy_output"
        }""")

        return default_settings


    def IsOutputStep(self):
        if self.output_control_is_time:
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
        # The name of the file refers to the step (as in vtk_output), independently of whether time or step is chosen for 'output_interval'.
        numpy.save(self.output_path / f"solution_{step}.npy", numpy.stack(aux_data_array, axis=1).reshape(-1,1))

        # Schedule the next output
        current = time if self.output_control_is_time else step
        while self.__GetPrettyFloat(self.next_output) <= self.__GetPrettyFloat(current):
            self.next_output += self.output_interval


    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number