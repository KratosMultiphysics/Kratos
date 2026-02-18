# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StatisticsApplication as ST
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

# Other imports
import time
from pathlib import Path

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RmsOutputProcess(model, settings["Parameters"])

class RmsOutputProcess(KM.Process):

    @staticmethod
    def GetDefaultParameters():
        """Return the default parameters."""
        return KM.Parameters("""
        {
            "model_part_name"         : "model_part",
            "filename"                : "output_file",
            "time_frequency"          : 1.0,
            "analysis_label"          : "label",
            "variables"               : [],
            "nonhistorical_variables" : []
        }
        """)


    def __init__(self, model, settings):
        """Constructor of RmsOutputProcess."""
        super().__init__()

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variables = GenerateVariableListFromInput(self.settings["variables"])
        self.nonhistorical_variables = GenerateVariableListFromInput(self.settings["nonhistorical_variables"])

        self.time_frequency = self.settings["time_frequency"].GetDouble()
        self.next_output = self.model_part.ProcessInfo[KM.TIME]

        self._InitializeOutputFile()


    def ExecuteInitialize(self):
        """Initialize the non historical variables and the measuring of computational time."""
        for variable in self.nonhistorical_variables:
            KM.VariableUtils().SetNonHistoricalVariableToZero(variable, self.model_part.Nodes)
        self.start_time = time.time()


    def IsOutputStep(self):
        """Check if the current time step is near enough to the specified printing times."""
        time = self.model_part.ProcessInfo[KM.TIME]
        return time >= self.next_output


    def PrintOutput(self):
        """Write the values into the file."""
        self._WriteAverageError()
        self.next_output += self.time_frequency


    def Check(self):
        """Check the correctness of the input parameters."""
        for variable in self.variables + self.nonhistorical_variables:
            if not isinstance(variable, KM.DoubleVariable):
                raise Exception("This process is expecting only double or component variables")


    def _InitializeOutputFile(self):
        file_path = Path(self.settings["filename"].GetString()).with_suffix('.dat')
        file_path.touch(exist_ok=True)
        header = self._GetFileHeader()
        if file_path.stat().st_size == 0:
            with open(file_path, 'w') as file:
                file.write(header)
        else:
            existing_header = ''
            with open(file_path, 'r') as file:
                for _ in range(2):
                    existing_header += file.readline()
            if existing_header != header:
                msg = self.__class__.__name__ + ": "
                msg += "The specified fields mismatch\n"
                msg += "Existing header:\n"
                msg += existing_header
                msg += "Specified header:\n"
                msg += header
                raise Exception(msg)


    def _GetFileHeader(self):
        header = "# RMS for model part '{}'\n".format(self.model_part.Name)
        header += "label datetime num_nodes num_elems elem_size time_step time computational_time"
        for variable in self.variables + self.nonhistorical_variables:
            header += ' ' + variable.Name()
        header += '\n'
        return header


    def _WriteAverageError(self):
        data  = self.settings["analysis_label"].GetString() + ' '
        data += str(time.time()) + ' '
        data += str(self.model_part.NumberOfNodes()) + ' '
        data += str(self.model_part.NumberOfElements()) + ' '
        data += str(self.model_part.Elements.__iter__().__next__().GetGeometry().Length()) + ' '
        data += str(self.model_part.ProcessInfo[KM.DELTA_TIME]) + ' '
        data += str(self.model_part.ProcessInfo[KM.TIME]) + ' '
        data += str(time.time() - self.start_time)

        for variable in self.variables:
            value = ST.SpatialMethods.Historical.NormMethods.RootMeanSquare(self.model_part, variable)
            data += ' {}'.format(value)
        for variable in self.nonhistorical_variables:
            value = ST.SpatialMethods.NonHistorical.Nodes.NormMethods.RootMeanSquare(self.model_part, variable)
            data += ' {}'.format(value)
        data += '\n'

        file_path = Path(self.settings["file_name"].GetString()).with_suffix('.dat')
        with open(file_path, 'a') as file:
            file.write(data)
