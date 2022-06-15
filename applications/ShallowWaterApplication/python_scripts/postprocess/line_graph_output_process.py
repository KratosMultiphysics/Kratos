import KratosMultiphysics as KM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput, DeleteFileIfExisting
from KratosMultiphysics.point_output_process import Interpolate
from pathlib import Path
from numpy import linspace


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return LineGraphOutputProcess(model, settings["Parameters"])


class LineGraphOutputProcess(KM.OutputProcess):
    """This process writes results along a line to generate a graph.
    Every output step, an output file will be generated containing the graph data.
    At the moment, MPI is not supported.
    """

    def GetDefaultParameters(self):
        return KM.Parameters("""
        {
            "help"                    : "This process writes results from a geometrical object (line) in the model to a file. It first searches the entities containing the requested output location and then interpolates the requested variable(s). The output can be requested for elements, conditions and nodes. For nodes no geometrical interpolation is performed, the exact coordinates have to be specified.",
            "model_part_name"         : "",
            "entity_type"             : "element",
            "interval"                : [0.0,"End"],
            "start_point"             : [0, 0, 0],
            "end_point"               : [0, 0, 0],
            "sampling_points"         : 100,
            "output_variables"        : [],
            "nonhistorical_variables" : [],
            "search_configuration"    : "initial",
            "search_tolerance"        : 1e-6,
            "print_format"            : "{:.6f}",
            "time_format"             : "{:.3f}",
            "output_file_settings"    : {
                "file_name"               : "<model_part>",
                "output_path"             : ""
            },
            "output_control_settings" : {
                "output_control_type"     : "time",
                "time_frequency"          : 1.0
            }
        }""")

    def __init__(self, model, settings):
        """Constructor of the class."""
        super().__init__()

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model_part = model[settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(settings)

        # Retrieving the positions defining the line entity
        start_point_position = settings["start_point"].GetVector()
        if start_point_position.Size() != 3:
            raise Exception('The start point position has to be provided with 3 coordinates!')
        end_point_position = settings["end_point"].GetVector()
        if end_point_position.Size() != 3:
            raise Exception('The end point position has to be provided with 3 coordinates!')

        # Get the number of points defining the line entity
        number_of_sampling_points = settings["sampling_points"].GetInt()
        if number_of_sampling_points <= 2:
            raise Exception('The number of sampling points has to be larger than 2!')
        parametrized_distances = linspace(0, 1, number_of_sampling_points)
        increment = end_point_position - start_point_position
        self.positions = [KM.Point(start_point_position + float(d)*increment) for d in parametrized_distances]

        # Check the entity type
        self.entity_type = settings["entity_type"].GetString()
        if not (self.entity_type == "node" or self.entity_type == "element" or self.entity_type == "condition"):
            raise Exception("Invalid 'entity_type' : {} (Expecting 'node', 'element' or 'condition')".format(self.entity_type))

        # Retrieve the variables list
        self.variables = self._GenerateVariablesList(settings["output_variables"], historical_value=True)
        self.nonhistorical_variables = self._GenerateVariablesList(settings["nonhistorical_variables"], historical_value=False)

        # Search settings
        if settings["search_configuration"].GetString() == "initial":
            self.search_configuration = KM.Configuration.Initial
        elif settings["search_configuration"].GetString() == "current":
            self.search_configuration = KM.Configuration.Current
        else:
            raise Exception("Invalid configuration: {} (Expecting 'initial' or 'current')".format(self.search_configuration))
        self.search_tolerance = settings["search_tolerance"].GetDouble()

        # Printing settings
        self.print_format = settings["print_format"].GetString()
        self.time_format = settings["time_format"].GetString()
        self.file_settings = settings["output_file_settings"].Clone()

        # Initialize output control
        self.output_control = OutputControlFactory(self.model_part, settings["output_control_settings"])


    def Check(self):
        """Check the file settings."""

        # Generate a dummy file to validate the parameters
        file = TimeBasedAsciiFileWriterUtility(self.model_part, self.file_settings, "").file
        file.close()
        DeleteFileIfExisting(file.name)


    def ExecuteBeforeSolutionLoop(self):
        """Search the points and delete the existing files after the current time."""

        # Get the file base name and check if there is a replacement
        self.file_name = self.file_settings["file_name"].GetString()
        self.file_name = self.file_name.replace("<model_part>", self.model_part.Name)

        # Delete the previous files
        time = self.model_part.ProcessInfo[KM.TIME]
        self._DeleteExistingFiles(time)

        # Perform the search
        self._SearchPoints()


    def IsOutputStep(self):
        """Return if the current step is an output step."""

        time = self.model_part.ProcessInfo[KM.TIME]
        return self.interval.IsInInterval(time) and self.output_control.IsOutputStep()


    def PrintOutput(self):
        """The output file is created, filled and closed.
        There will be one file for each printing step with the time as label.
        """

        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        dummy_extension = ".z" #NOTE: the dummy extension will be replaced by the file utility. It is used to keep the decimals.
        self.file_settings["file_name"].SetString(self.file_name + '_' + self.time_format.format(time) + dummy_extension)
        file = TimeBasedAsciiFileWriterUtility(self.model_part, self.file_settings, self._GetHeader()).file
        for point, entity, area_coords in zip(self.found_positions, self.entities, self.area_coords):
            file.write(self._GetPointData(point, entity, area_coords))
        file.close()


    def _GenerateVariablesList(self, parameters, historical_value):
        all_variables_list = GenerateVariableListFromInput(parameters)
        variables = []
        # Validate the types of variables
        for var in all_variables_list:
            if historical_value:
                if not self.model_part.HasNodalSolutionStepVariable(var):
                    raise Exception("ModelPart '{}' does not have {} as SolutionStepVariable".format(self.model_part.Name, var.Name()))
            if isinstance(var, KM.DoubleVariable):
                variables.append(var)
            elif isinstance(var, KM.Array1DVariable3):
                variables.append(KM.KratosGlobals.GetVariable(var.Name() + "_X"))
                variables.append(KM.KratosGlobals.GetVariable(var.Name() + "_Y"))
                variables.append(KM.KratosGlobals.GetVariable(var.Name() + "_Z"))
            else:
                raise Exception("The variable {} is not valid. It can only be double, component or array_3d".format(var.Name()))
        return variables


    def _SearchPoints(self):
        self.entities = []
        self.area_coords = []
        self.found_positions = []
        if self.entity_type == "node":
            for point in self.positions:
                found_id = KM.BruteForcePointLocator(self.model_part).FindNode(point, self.search_configuration, self.search_tolerance)
                if found_id > -1:
                    self.entities.append(self.model_part.Nodes[found_id])
                    self.area_coords.append("dummy") # needed for looping later
                    self.found_positions.append(point)
        elif self.entity_type == "element":
            for point in self.positions:
                self.sf_values = KM.Vector()
                found_id = KM.BruteForcePointLocator(self.model_part).FindElement(point, self.sf_values, self.search_configuration, self.search_tolerance)
                if found_id > -1:
                    self.entities.append(self.model_part.Elements[found_id])
                    self.area_coords.append(self.sf_values)
                    self.found_positions.append(point)
        elif self.entity_type == "condition":
            for point in self.positions:
                self.sf_values = KM.Vector()
                found_id = KM.BruteForcePointLocator(self.model_part).FindCondition(point, self.sf_values, self.search_configuration, self.search_tolerance)
                if found_id > -1:
                    self.entities.append(self.model_part.Conditions[found_id])
                    self.area_coords.append(self.sf_values)
                    self.found_positions.append(point)


    def _GetHeader(self):
        if len(self.found_positions) > 1:
            start = list(self.found_positions[0])
            end = list(self.found_positions[-1])
        else:
            start = "'NOT FOUND'"
            end = "'NOT FOUND'"
        time = self.model_part.ProcessInfo[KM.TIME]
        header = "# Results for '{}s' over line {}-{} at time {}\n#".format(self.entity_type, start, end, time)
        coordinates = ["X", "Y", "Z"]
        for c in coordinates:
            header += " " + c
        for var in self.variables:
            header += " " + var.Name()
        for var in self.nonhistorical_variables:
            header += " " + var.Name()
        return header + "\n"


    def _GetPointData(self, node, entity, area_coords):
        data = self.print_format.format(node.X)
        data += " " + self.print_format.format(node.Y)
        data += " " + self.print_format.format(node.Z)
        for var in self.variables:
            data += " " + self.print_format.format(Interpolate(var, entity, area_coords, historical_value=True))
        for var in self.nonhistorical_variables:
            data += " " + self.print_format.format(Interpolate(var, entity, area_coords, historical_value=False))
        return data + "\n"


    def _DeleteExistingFiles(self, time):
        output_path = self.file_settings["output_path"].GetString()
        file_extension = self.file_settings["file_extension"].GetString()
        for dir in Path(output_path).glob(self.file_name + "*." + file_extension):
            file_time = dir.stem.split("-")[-1]
            try:
                file_time = float(file_time)
                if file_time >= time:
                    DeleteFileIfExisting(dir)
            except ValueError:
                pass




# TODO: move this to a common place

def OutputControlFactory(model_part, parameters):
    if not parameters.Has("output_control_type"):
        raise Exception ("OutputControlFactory. There is no 'output_control_settings' key in the parameters.")
    output_control_type = parameters["output_control_type"].GetString()

    output_control_types = {
        "time" : TimeOutputControl,
        "step" : StepOutputControl
    }
    module = output_control_types[output_control_type]
    return module(model_part, parameters)

class OutputControl():

    @staticmethod
    def IsOutputStep():
        return True

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters()

class TimeOutputControl(OutputControl):

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "output_control_type"     : "time",
            "time_frequency"          : 1.0
        }""")

    def __init__(self, model_part, parameters):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.time_frequency = parameters["time_frequency"].GetDouble()
        self.model_part = model_part
        self.next_output = 0.0

    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KM.TIME]
        if time >= self.next_output:
            while time >= self.next_output:
                self.next_output += self.time_frequency
            return True
        return False

class StepOutputControl(OutputControl):

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "output_control_type"     : "step",
            "step_frequency"          : 1
        }""")

    def __init__(self, model_part, parameters):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.step_frequency = parameters["step_frequency"].GetDouble()
        self.model_part = model_part
        self.next_output = 0

    def IsOutputStep(self):
        step = self.model_part.ProcessInfo[KM.STEP]
        if step >= self.next_output:
            while step >= self.next_output:
                self.next_output += self.step_frequency
            return True
        return False
