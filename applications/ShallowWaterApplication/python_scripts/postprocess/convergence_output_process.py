# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

# Other imports
import time
from pathlib import Path

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ConvergenceOutputProcess(model, settings["Parameters"])

class ConvergenceOutputProcess(KM.Process):

    def __init__(self, model, settings):
        super().__init__()

        default_settings = KM.Parameters("""{
                "model_part_name"       : "model_part",
                "file_name"             : "output_file",
                "printing_times"        : [],
                "analysis_label"        : "label",
                "convergence_variables" : [],
                "low_corner"            : [],
                "high_corner"           : []
            }""")

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variables = GenerateVariableListFromInput(self.settings["convergence_variables"])

        # Initialize output control variables
        self.printing_times = self.settings["printing_times"].GetVector()
        self.is_printed = [False] * len(self.printing_times)

        if self.settings["low_corner"].GetVector().Size() == 0:
            self.integrate_over_all_the_domain = True
        else:
            self.integrate_over_all_the_domain = False

        self.file = self._InitializeOutputFile()
        self.file.flush() # Since this process will reuse the existing files, we need to ensure the file construction is finished. Otherwise problems may arise if the first analysis didn't finish.


    def ExecuteInitialize(self):
        """Initialize the non historical variables and the measuring of computational time."""
        for variable in self.variables:
            KM.VariableUtils().SetNonHistoricalVariableToZero(variable, self.model_part.Nodes)
        self.start_time = time.time()


    def IsOutputStep(self):
        """Check if the current time step is near enough to the specified printing times."""
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        for i in range(len(self.printing_times)):
            if time >= self.printing_times[i] and not self.is_printed[i]:
                self.is_printed[i] = True
                return True
        return False


    def PrintOutput(self):
        """Write the values into the file."""
        self._WriteAverageError()


    def ExecuteFinalize(self):
        """Close the file."""
        self.file.close()


    def Check(self):
        """Check the correctness of the input parameters."""
        for variable in self.variables:
            if not isinstance(variable, KM.DoubleVariable):
                raise Exception("This process is expecting only double or component variables")
        
        low_corner = self.settings["low_corner"].GetVector()
        high_corner = self.settings["high_corner"].GetVector()
        if not low_corner.Size() == high_corner.Size():
            raise Exception("The low and high corners does not have the same dimension")

        if low_corner.Size() == 0:
            pass
        elif low_corner.Size() == 2:
            self.settings["low_corner"].Append(0.0)
            self.settings["high_corner"].Append(0.0)
        elif low_corner.Size() == 3:
            pass
        else:
            raise Exception("The corners must be specified with 2 or 3 coordinates")


    def _InitializeOutputFile(self):
        path = Path(self.settings["file_name"].GetString()).with_suffix('.dat')
        file = open(path, "a+")
        if path.stat().st_size == 0:
            header = self._GetHeader()
            file.write(header)
        return file


    def _GetHeader(self):
        header = "# RMS for model part '{}' ".format(self.model_part.Name)
        if self.integrate_over_all_the_domain:
            header += "over all the domain\n"
        else:
            low_corner = KM.Point(self.settings["low_corner"].GetVector())
            high_corner = KM.Point(self.settings["high_corner"].GetVector())
            header += "over rectangle {} x {}\n".format(list(low_corner), list(high_corner))

        header += "#label \t num_nodes \t num_elems \t time_step \t time \t computational_time"
        for variable in self.variables:
            header += '\t' + variable.Name()
        header += '\n'
        return header


    def _WriteAverageError(self):
        data  = self.settings["analysis_label"].GetString() + '\t'
        data += str(self.model_part.NumberOfNodes()) + '\t'
        data += str(self.model_part.NumberOfElements()) + '\t'
        data += str(self.model_part.ProcessInfo[KM.DELTA_TIME]) + '\t'
        data += str(self.model_part.ProcessInfo[KM.TIME]) + '\t'
        data += str(time.time() - self.start_time)

        if not self.integrate_over_all_the_domain:
            low_corner = KM.Point(self.settings["low_corner"].GetVector())
            high_corner = KM.Point(self.settings["high_corner"].GetVector())

        for variable in self.variables:
            if self.integrate_over_all_the_domain:
                value = SW.ShallowWaterUtilities().ComputeL2NormNonHistorical(self.model_part, variable)
            else:
                value = SW.ShallowWaterUtilities().ComputeL2NormNonHistorical(self.model_part, variable, low_corner, high_corner)
            data += '\t{}'.format(value)
        data += '\n'

        self.file.write(data)
        self.file.flush()
