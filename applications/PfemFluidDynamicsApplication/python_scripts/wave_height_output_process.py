import KratosMultiphysics as KM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveHeightOutputProcess(model, settings["Parameters"])

class WaveHeightOutputProcess(KM.OutputProcess):
    """WaveHeightOutputProcess

    This process records the wave height at several points.
    The direction used to calculate the water weight is defined as the opposite of the gravity direction.
    If no node is found, Nan will be printed.
    Possible specifications of the Parameters:
     - coordinates: it can be a single coordinate or a list of coordinates for each gauge.
     - file_names:  it can be a single string or a list of string. If there is a single string and
                    multiple coordinates, the string will be repeated for all the coordinates.
                    Some replacements can be added, such as 'gauge_<X>'.
     - output_path: same as file_names.
    """

    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""{
            "model_part_name"        : "",
            "coordinates"            : [[0.0, 0.0, 0.0]],
            "mean_water_level"       : 0.0,
            "relative_search_radius" : 1.0,
            "search_tolerance"       : 1e-6,
            "file_names"             : [""],
            "output_path"            : [""],
            "time_between_outputs"   : 0.01
        }""")
        if self.settings.Has("file_names"):
            if self.settings["file_names"].IsString():
                default_parameters["file_names"].SetString("")
        if self.settings.Has("output_path"):
            if self.settings["output_path"].IsString():
                default_parameters["output_path"].SetString("")
        return default_parameters

    def __init__(self, model, settings):
        """The constructor of the WaveHeightOutputProcess"""

        KM.OutputProcess.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())

        self.coordinates_list = self._GetCoordinatesList()
        self.file_names_list = self._GetNamesList(self.settings["file_names"])
        self.output_path_list = self._GetNamesList(self.settings["output_path"])

        self.next_output = self.model_part.ProcessInfo[KM.TIME]

    def Check(self):
        """Check all the variables have the right size"""
        for coordinate in self.coordinates_list:
            if not len(coordinate) == 3:
                raise Exception("WaveHeightOutputProcess. The coordinates must be given with a 3-dimensional array")

        if not len(self.file_names_list) == len(self.coordinates_list):
            raise Exception("WaveHeightOutputProcess. The number of coordinates must coincide with the number of filenames")

        if not len(self.output_path_list) == len(self.coordinates_list):
            raise Exception("WaveHeightOutputProcess. The number of coordinates must coincide with the number of output paths")

    def ExecuteBeforeSolutionLoop(self):
        """Initialize the files and the utility to calculate the water height"""
        self.files = []
        for file_name, output_path in zip(self.file_names_list, self.output_path_list):
            file_settings = KM.Parameters()
            file_settings.AddString("file_name", file_name)
            file_settings.AddString("output_path", output_path)
            header = "#Time \tHeight\n"
            self.files.append(TimeBasedAsciiFileWriterUtility(self.model_part, file_settings, header).file)
        utility_settings = KM.Parameters()
        utility_settings.AddValue("mean_water_level", self.settings["mean_water_level"])
        utility_settings.AddValue("search_tolerance", self.settings["search_tolerance"])
        utility_settings.AddValue("relative_search_radius", self.settings["relative_search_radius"])
        self.wave_height_utility = PFEM.CalculateWaveHeightUtility(self.model_part, utility_settings)

    def IsOutputStep(self):
        """The output control is time based"""
        time = self.model_part.ProcessInfo[KM.TIME]
        return time >= self.next_output

    def PrintOutput(self):
        """Print the wave height corresponding to each gauge and schedule the next output"""
        time = self.model_part.ProcessInfo[KM.TIME]
        for file, coordinates in zip(self.files, self.coordinates_list):
            height = self.wave_height_utility.Calculate(coordinates)
            value = "{}\t{}\n".format(time, height)
            file.write(value)
            file.flush()

        self.next_output += self.settings["time_between_outputs"].GetDouble()

    def ExecuteFinalize(self):
        """Close all the files"""
        for file in self.files:
            file.close()

    def _GetCoordinatesList(self):
        coordinates_list = []
        if self.settings["coordinates"].IsVector(): # There is a single coordinate
            coordinates_list.append(self.settings["coordinates"].GetVector())
        else:
            for coordinate in self.settings["coordinates"]: # There is a list of coordinates
                coordinates_list.append(coordinate.GetVector())
        return coordinates_list

    def _GetNamesList(self, settings):
        if settings.IsString():
            names_list = []
            for coordinate in self.coordinates_list:
                name = settings.GetString()
                name = name.replace("<x>", str(coordinate[0]))
                name = name.replace("<X>", str(coordinate[0]))
                name = name.replace("<y>", str(coordinate[1]))
                name = name.replace("<Y>", str(coordinate[1]))
                name = name.replace("<z>", str(coordinate[2]))
                name = name.replace("<Z>", str(coordinate[2]))
                names_list.append(name)
            return names_list
        else:
            return settings.GetStringArray()
