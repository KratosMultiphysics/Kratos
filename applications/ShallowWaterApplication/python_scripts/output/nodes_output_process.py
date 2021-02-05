import KratosMultiphysics as KM
import numpy as np
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return NodesOutputProcess(model, settings["Parameters"])

class NodesOutputProcess(KM.Process):
    """This process writes results from the nodes which are near to a line.
    If the distance from a node to the line is less than a tolerance, it is added
    to the output.
    The line is extended from the start and ending point.
    """
    def __init__(self, model, settings):
        """Constructor of the class."""
        super().__init__()

        default_settings = KM.Parameters("""
            {
                "help"                    : "This process writes results from the nodes which are near a line.",
                "model_part_name"         : "model_part",
                "file_name"               : "output_file",
                "output_path"             : "",
                "start_point"             : [],
                "end_point"               : [],
                "output_variables"        : [],
                "nonhistorical_variables" : [],
                "printing_times"          : [],
                "write_buffer_size"       : -1,
                "relative_tol_to_line"    : 0.25
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]

        # Retrieving the positions defining the line entity
        start_point_position = self.settings["start_point"].GetVector()
        if start_point_position.Size() != 3:
            raise Exception('The start point position has to be provided with 3 coordinates!')
        end_point_position = self.settings["end_point"].GetVector()
        if end_point_position.Size() != 3:
            raise Exception('The end point position has to be provided with 3 coordinates!')
        self.origin = start_point_position
        self.direction = end_point_position - self.origin

        # The output file settings
        self.out_file_params = KM.Parameters()
        self.out_file_params.AddEmptyValue("file_name")
        self.out_file_params.AddValue("output_path", self.settings["output_path"])
        self.out_file_params.AddValue("write_buffer_size", self.settings["write_buffer_size"])

        # The tolerance relative tho the average element size
        self.rel_tolerance = self.settings["relative_tol_to_line"].GetDouble()

    def ExecuteInitialize(self):
        """The model part related variables are initialized:
        The list of nodes and the list of variables are created.
        """
        self.nodes = []
        tolerance = self._GetTolerance()
        for node in self.model_part.Nodes:
            if self._DistanceToLine(node) < tolerance:
                self.nodes.append(node)

        self.variables, self.names = self._GenerateVariablesList(self.settings["output_variables"])
        self.nonhist_variables, self.nonhist_names = self._GenerateVariablesList(self.settings["nonhistorical_variables"], False)

    def Check(self):
        """This check verifies if the two specified points are found inside the domain"""
        start_point_position = self.settings["start_point"].GetVector()
        end_point_position = self.settings["end_point"].GetVector()
        start_point = KM.Point(start_point_position)
        end_point = KM.Point(end_point_position)
        locator = KM.BruteForcePointLocator(self.model_part)
        if locator.FindNode(start_point, KM.Configuration.Initial, self._GetTolerance()) < 0:
            raise Exception('The start point was not found in the domain. Please, check the geometry or the relative tolerance')
        if locator.FindNode(end_point, KM.Configuration.Initial, self._GetTolerance()) < 0:
            raise Exception('The end point was not found in the domain. Please, check the geometry or the relative tolerance')

    def IsOutputStep(self):
        """This method checks if the current time step is
        near enough to the specified printing times.
        """
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        delta_time = self.model_part.ProcessInfo.GetValue(KM.DELTA_TIME)
        for printing_time in self.settings["printing_times"].GetVector():
            if 2 * abs(time - printing_time) < delta_time:
                return True
        return False

    def PrintOutput(self):
        """The output file is created, filled and closed. If several output
        timesteps are specified, there will be one file for each timestep.
        """
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        file_name = self.settings["file_name"].GetString() + "_{:.4f}.dat".format(time)
        self.out_file_params["file_name"].SetString(file_name)
        file = TimeBasedAsciiFileWriterUtility(self.model_part, self.out_file_params, self._GetHeader()).file
        for node in self.nodes:
            file.write(self._GetData(node, self._DistanceToOrigin(node)))
        file.close()

    def _GetTolerance(self):
        if not hasattr(self, 'tolerance'):
            elem_size = 0
            for elem in self.model_part.Elements:
                elem_size += elem.GetGeometry().Length()
            elem_size /= self.model_part.NumberOfElements()
            self.tolerance = self.rel_tolerance * elem_size
        return self.tolerance

    def _DistanceToLine(self, point):
        return np.linalg.norm(np.cross(self.direction, point - self.origin)) / np.linalg.norm(self.direction)

    def _DistanceToOrigin(self, point):
        return np.linalg.norm(point - self.origin)

    def _GenerateVariablesList(self, parameters, historical=True):
        variables_names = parameters.GetStringArray()
        names = []
        variables = []
        for name in variables_names:
            var = KM.KratosGlobals.GetVariable(name)
            if historical:
                self._CheckIsHistoricalVariable(var)
            self._AppendVariable(var, variables, names)
        return variables, names

    def _CheckIsHistoricalVariable(self, input_variable):
        if not self.model_part.HasNodalSolutionStepVariable(input_variable):
            err_msg  = 'ModelPart "' + self.model_part.Name + '" does not have'
            err_msg += ' "' + input_variable.Name + '" as SolutionStepVariable!'
            raise Exception(err_msg)

    @staticmethod
    def _AppendVariable(input_variable, variables, names):
        name = input_variable.Name()
        if isinstance(input_variable, KM.DoubleVariable):
            names.append(name)
            variables.append(input_variable)
        elif isinstance(input_variable, KM.Array1DVariable3):
            names.append(name + "_X")
            names.append(name + "_Y")
            names.append(name + "_Z")
            variables.append(KM.KratosGlobals.GetVariable(name + "_X"))
            variables.append(KM.KratosGlobals.GetVariable(name + "_Y"))
            variables.append(KM.KratosGlobals.GetVariable(name + "_Z"))
        else:
            err_msg  = 'Type of variable "' + name + '" is not valid\n'
            err_msg += 'It can only be double, component or array3d!'
            raise Exception(err_msg)

    def _GetHeader(self):
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        header  = "# Nodes along line:\n"
        header += "#    origin = ({:.4f}, {:.4f}, {:.4f})\n".format(*self.origin)
        header += "#    direction = ({:.4f}, {:.4f}, {:.4f})\n".format(*self.direction)
        header += "# Time = {:.4f}\n".format(time)
        header += "#\n"
        header += "#position "
        for name in self.names:
            header += name + " "
        for name in self.nonhist_names:
            header += name + " "
        return header + "\n"

    def _GetData(self, node, position):
        line = str(position) + " "
        for var in self.variables:
            line += str(node.GetSolutionStepValue(var)) + " "
        for var in self.nonhist_variables:
            line += str(node.GetValue(var)) + " "
        return line + "\n"
