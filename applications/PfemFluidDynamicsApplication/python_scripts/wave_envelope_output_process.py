import KratosMultiphysics as KM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from numpy import linspace

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveHeightOutputProcess(model, settings["Parameters"])

class WaveHeightOutputProcess(KM.OutputProcess):
    """WaveHeightOutputProcess

    This process records the wave height along a path.
    If a sampling point is not found, Nan will be printed.
    Possible specifications of the Parameters:
     - output_file_settings: a parameters encapsulating the 'file_name', 'output_path' and
                             other settings according to 'TimeBasedAsciiFileWritterUtility'.
     - wave_calculation_settings: a parameters according to 'CalculateWaveHeightOutputProcess'
    """

    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""{
            "model_part_name"           : "",
            "start_point"               : [0, 0, 0],
            "end_point"                 : [0, 0, 0],
            "sampling_points"           : 100,
            "wave_calculation_settings" : {},
            "output_file_settings"      : {},
            "time_between_outputs"      : 1.0,
            "print_format"              : "{:.6f}"
        }""")
        return default_parameters


    def __init__(self, model, settings):
        """The constructor of the WaveHeightOutputProcess"""

        KM.OutputProcess.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())

        # Retrieving the positions defining the line entity
        start_point_position = settings["start_point"].GetVector()
        end_point_position = settings["end_point"].GetVector()

        # Get the number of points defining the line entity
        number_of_sampling_points = settings["sampling_points"].GetInt()
        if number_of_sampling_points <= 2:
            raise Exception('The number of sampling points has to be larger than 2!')
        parametrized_distances = linspace(0, 1, number_of_sampling_points)
        abs_increment = end_point_position - start_point_position
        self.positions = [KM.Point(start_point_position + float(d)*abs_increment) for d in parametrized_distances]

        # Schedule the next output
        self.next_output = self.model_part.ProcessInfo[KM.TIME]


    def ExecuteBeforeSolutionLoop(self):
        """Initialize the utility to calculate the water height"""
        wave_settings = self.settings["wave_calculation_settings"]
        self.wave_height_utility = PFEM.CalculateWaveHeightUtility(self.model_part, wave_settings)
        self.max_values = [-1e6] * len(self.positions)


    def ExecuteFinalizeSolutionStep(self):
        """Calculate the maximum wave height along the path."""
        for i, point in enumerate(self.positions):
            value = self.wave_height_utility.Calculate(point)
            self.max_values[i] = max(self.max_values[i], value)


    def IsOutputStep(self):
        """The output control is time based"""
        time = self.model_part.ProcessInfo[KM.TIME]
        return time >= self.next_output


    def PrintOutput(self):
        """Print the wave height corresponding to each gauge and schedule the next output.

        The previous output files are overwritten. If the simulation
        does not reach the end, an envelope will be kept.
        """
        file_settings = self.settings["output_file_settings"]
        header = "# Wave height along line {}-{}\n"
        header += "# X \t\t Y \t\t Z \t\t Height\n"
        header = header.format(list(self.positions[0]), list(self.positions[-1]))
        file = TimeBasedAsciiFileWriterUtility(self.model_part, file_settings, header).file

        print_format = self.settings["print_format"].GetString()
        row = print_format + "\t" + print_format + "\t" + print_format + "\t" + print_format + "\n"
        for point, value in zip(self.positions, self.max_values):
            value = row.format(point.X, point.Y, point.Z, value)
            file.write(value)
        file.close()

        # Schedule the next output
        self.next_output += self.settings["time_between_outputs"].GetDouble()
