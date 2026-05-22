# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
import math

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a model object")

    return MPMParticleHeightOutputProcess(model, settings["Parameters"])


class MPMParticleHeightOutputProcess(KratosMultiphysics.OutputProcess):
    """This process writes the height of material particles close to sensor lines."""

    def __init__(self, model, params):
        KratosMultiphysics.OutputProcess.__init__(self)

        default_settings = KratosMultiphysics.Parameters('''{
            "help"                            : "This process writes the maximum particle height measured along a direction for one or more sensor lines.",
            "model_part_name"                 : "",
            "background_grid_model_part_name" : "",
            "interval"                        : [0.0, 1e30],
            "sensor_positions"                : [],
            "measuring_direction"             : [0.0, 1.0, 0.0],
            "search_radius_factor"            : 0.5,
            "print_format"                    : ".8f",
            "output_file_settings"            : {}
        }''')

        self.model = model
        self.params = params

        if self.params.Has("interval"):
            if self.params["interval"][1].IsString():
                if self.params["interval"][1].GetString() == "End":
                    self.params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" + self.params["interval"].PrettyPrintJsonString())

        self.params.ValidateAndAssignDefaults(default_settings)

        self.model_part_name = self.params["model_part_name"].GetString()
        if self.model_part_name == "":
            raise Exception('No "model_part_name" was specified!')

        self.model_part = self.model[self.model_part_name]
        self.format = self.params["print_format"].GetString()
        self.output_file = None

    def ExecuteBeforeSolutionLoop(self):
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = self.params["interval"][0].GetDouble()
        self.interval[1] = self.params["interval"][1].GetDouble()

        self.sensor_positions = self._ReadSensorPositions()
        self.measuring_direction = self._NormalizeVector(self.params["measuring_direction"].GetVector(), "measuring_direction")
        self.search_radius = self._CalculateSearchRadius()

        file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])
        if not file_handler_params.Has("file_name"):
            output_file_name = self.model_part_name + "_particle_height"
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "File name not specified, using the default name: " + output_file_name)
            file_handler_params.AddString("file_name", output_file_name)

        file_header = self._GetFileHeader()
        self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_params, file_header).file

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if current_time >= self.interval[0] and current_time <= self.interval[1]:
            heights = self._CalculateHeights()
            self._WriteOutput(current_time, heights)

    def IsOutputStep(self):
        return False

    def PrintOutput(self):
        pass

    def ExecuteFinalize(self):
        if self.output_file:
            self.output_file.close()

    def _ReadSensorPositions(self):
        sensor_positions = []
        settings_sensor_positions = self.params["sensor_positions"]

        if settings_sensor_positions.size() == 0:
            raise Exception('At least one position has to be provided in "sensor_positions"!')

        for i in range(settings_sensor_positions.size()):
            position = settings_sensor_positions[i].GetVector()
            if position.Size() != 3:
                raise Exception('Each sensor position has to be provided with 3 coordinates!')
            sensor_positions.append(position)

        return sensor_positions

    def _CalculateSearchRadius(self):
        background_grid_model_part_name = self.params["background_grid_model_part_name"].GetString()
        if background_grid_model_part_name == "":
            raise Exception('A "background_grid_model_part_name" has to be provided to calculate the search radius!')

        background_grid_model_part = self.model[background_grid_model_part_name]
        mesh_size = self._EstimateMeshSize(background_grid_model_part)
        return self.params["search_radius_factor"].GetDouble() * mesh_size

    def _EstimateMeshSize(self, model_part):
        number_of_elements = 0
        total_element_size = 0.0

        for element in model_part.Elements:
            element_size = element.GetGeometry().Length()
            if element_size > 0.0:
                number_of_elements += 1
                total_element_size += element_size

        if number_of_elements == 0:
            raise Exception("The background grid model part does not have enough geometry to estimate the mesh size!")

        return total_element_size / number_of_elements

    def _CalculateHeights(self):
        heights = []
        for sensor_position in self.sensor_positions:
            heights.append(float("nan"))

        for mp in self.model_part.Elements:
            mp_coord = mp.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, self.model_part.ProcessInfo)[0]

            for i in range(len(self.sensor_positions)):
                distance_to_sensor = self._DistanceToSensorLine(mp_coord, self.sensor_positions[i], self.measuring_direction)
                if distance_to_sensor <= self.search_radius:
                    height = self._Dot(mp_coord - self.sensor_positions[i], self.measuring_direction)
                    if math.isnan(heights[i]) or height > heights[i]:
                        heights[i] = height

        return heights

    def _DistanceToSensorLine(self, point, sensor_position, direction):
        sensor_to_point = point - sensor_position
        height = self._Dot(sensor_to_point, direction)
        closest_point = sensor_position + height * direction
        return (point - closest_point).norm_2()

    def _WriteOutput(self, time, heights):
        output = str(time)
        for height in heights:
            output += " " + format(height, self.format)
        self.output_file.write(output + "\n")
        self.output_file.flush()

    def _GetFileHeader(self):
        header  = "# Particle height measured on model part '" + self.model_part_name + "'\n"
        header += "# Search radius: " + str(self.search_radius) + "\n"
        header += "# Measuring direction: " + str(self.measuring_direction) + "\n"
        header += "# time"
        for i in range(len(self.sensor_positions)):
            header += " height_" + str(i + 1)
        header += "\n"
        return header

    def _NormalizeVector(self, vector, vector_name):
        if vector.Size() != 3:
            raise Exception('The "' + vector_name + '" has to be provided with 3 coordinates!')

        norm = vector.norm_2()
        if norm == 0.0:
            raise Exception('The "' + vector_name + '" cannot be a zero vector!')

        return vector / norm

    def _Dot(self, a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
