import KratosMultiphysics as KM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.PfemFluidDynamicsApplication.wave_height_output_process import WaveHeightOutputProcess

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ClosestWaveHeightOutputProcess(model, settings["Parameters"])

class ClosestWaveHeightOutputProcess(WaveHeightOutputProcess):
    """ClosestWaveHeightOutputProcess

    This process records the wave height at several points.
    The direction used to calculate the water weight is defined as the opposite of the gravity direction.
    If no node is found, Nan will be printed.
    Possible specifications of the Parameters:
     - coordinates: it can be a single coordinate or a list of coordinates for each gauge.
     - mesh_size: list with the mesh size around each gauge. If not defined, then the average size of all elements is used
     - output_file_settings: a parameters encapsulating the 'file_name', 'output_path' and
                             other settings according to 'TimeBasedAsciiFileWritterUtility'.
                             Some replacements can be specified to 'file_name', e.g.:
                              - 'file_name' : 'gauge_<x>' or
                              - 'file_name' : 'gauge_<Y>' or
                              - 'file_name' : 'gauge_<i>'
     - wave_calculation_settings: a parameters according to 'CalculateClosestWaveHeightOutputProcess'
                              - 'mean_water_level'
    """

    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""{
            "model_part_name"           : "",
            "coordinates"               : [[0.0, 0.0, 0.0]],
            "wave_calculation_settings" : {},
            "output_file_settings"      : {},
            "time_between_outputs"      : 0.01,
            "print_format"              : "{:.6f}"
        }""")
        return default_parameters

    def __init__(self, model, settings):
        """The constructor of the ClosestWaveHeightOutputProcess"""

        KM.OutputProcess.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())
        self.coordinates_list = self._GetCoordinatesList(self.settings["coordinates"])
        self.next_output = self.model_part.ProcessInfo[KM.TIME]

    def Check(self):
        """Check all the variables have the right size"""
        for coordinate in self.coordinates_list:
            if not len(coordinate) == 3:
                raise Exception("ClosestWaveHeightOutputProcess. The coordinates must be given with a 3-dimensional array")

        if len(self.coordinates_list) > 1:
            # Check 'file_name' exists and there is a possible replacement for the multiple output locations
            file_name = self.settings["output_file_settings"]["file_name"].GetString()
            if not '<' in file_name or not '>' in file_name:
                file_name += '_<i>'
                self.settings["output_file_settings"]["file_name"].SetString(file_name)

    def ExecuteBeforeSolutionLoop(self):
        """Initialize the files and the utility to calculate the water height"""
        self.files = []
        for i, coordinate in enumerate(self.coordinates_list, start=1):
            file_settings = self.settings["output_file_settings"].Clone()
            self._ExecuteReplacement(i, coordinate, file_settings["file_name"])
            header = "# Wave height at coordinate {}\n".format(list(coordinate))
            header += "#Time \t\tHeight\n"
            self.files.append(TimeBasedAsciiFileWriterUtility(self.model_part, file_settings, header).file)

        ## this gives the height of the closest node to the gauge
        self.wave_height_utility = PFEM.CalculateClosestWaveHeightUtility(self.model_part, self.settings["wave_calculation_settings"])

    def PrintOutput(self):
        """Print the wave height corresponding to each gauge and schedule the next output"""
        print_format = self.settings["print_format"].GetString()
        row = print_format + "\t" + print_format + "\n"
        time = self.model_part.ProcessInfo[KM.TIME]
        for file, coordinates in zip(self.files, self.coordinates_list):
            height = self.wave_height_utility.Calculate(coordinates)
            value = row.format(time, height)
            file.write(value)
            file.flush()

        self.next_output += self.settings["time_between_outputs"].GetDouble()

