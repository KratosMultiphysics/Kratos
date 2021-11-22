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
    """

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "model_part_name"      : ""
            "coordinates_list"     : [[0.0, 0.0, 0.0]],
            "mean_water_level"     : 0.0,
            "search_tolerance"     : 1.0,
            "file_names_list"      : [""],
            "output_path"          : "",
            "time_between_outputs" : 0.01
        }""")

    def __init__(self, model, settings):
        """The constructor of the WaveHeightOutputProcess"""

        KM.Process.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())

        self.coordinates_list = []
        for coordinate in self.settings["coordinates_list"]:
            self.coordinates_list.append(coordinate.GetVector())
        self.file_names_list = self.settings["file_names_list"]

    def Check(self):
        for coordinate in self.coordinates_list:
            if not coordinate.size() == 3:
                raise Exception("WaveHeightOutputProcess. The coordinates must be given with a 3-dimensional array")

        if not self.file_names_list.size() == self.coordinates_list.size():
            raise Exception("WaveHeightOutputProcess. The number of coordinates must coincide with the number of filenames")

    def ExecuteBeforeSolutionLoop(self):
        self.files = []
        for file_name in self.file_names_list:
            file_settings = KM.Parameters()
            file_settings.AddString("file_name", file_name)
            file_settings.AddValue("output_path", self.settings["output_path"])
            header = "#Time \tHeight\n"
            self.files.append(TimeBasedAsciiFileWriterUtility(self.model_part, file_settings, header).file)
        utility_settings = KM.Parameters()
        utility_settings.AddValue("mean_water_level", self.settings["mean_water_level"])
        utility_settings.AddValue("search_tolerance", self.settings["search_tolerance"])
        self.wave_height_utility = PFEM.CalculateWaveHeightUtility(self.model_part, utility_settings)

    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KM.TIME]
        return time >= self.next_output
    
    def PrintOutput(self):
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        for file, coordinates in zip(self.files, self.coordinates_list):
            height = self.wave_height_utility.Calculate(coordinates)
            value = "{}\t{}\n".format(time, height)
            file.write(value)
            file.flush()
        
        self.next_output += self.settings["time_between_outputs"].GetDouble()
    
    def ExecuteFinalize(self):
        for file in self.files:
            file.close()
