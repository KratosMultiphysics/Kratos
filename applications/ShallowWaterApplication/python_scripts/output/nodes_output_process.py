import KratosMultiphysics as KM
import numpy as np
import os, glob
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
                "relative_tol_to_line"    : 0.1,
                "use_mesh_nodes"          : true,
                "number_of_nodes"         : 200
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]

        # Retrieving the positions defining the line entity
        start_point_position = self.settings["start_point"].GetVector()
        if start_point_position.Size() == 2:
            self.settings["start_point"].Append(0.0)
        elif start_point_position.Size() != 3:
            raise Exception('The start point position has to be provided with 2 or 3 coordinates!')

        end_point_position = self.settings["end_point"].GetVector()
        if end_point_position.Size() == 2:
            self.settings["end_point"].Append(0.0)
        elif end_point_position.Size() != 3:
            raise Exception('The end point position has to be provided with 2 or 3 coordinates!')

        self.origin = self.settings["start_point"].GetVector()
        self.direction = self.settings["end_point"].GetVector() - self.origin

        # The output file settings
        self.out_file_params = KM.Parameters()
        self.out_file_params.AddEmptyValue("file_name")
        self.out_file_params.AddValue("output_path", self.settings["output_path"])
        self.out_file_params.AddValue("write_buffer_size", self.settings["write_buffer_size"])

        # The tolerance relative tho the average element size
        self.rel_tolerance = self.settings["relative_tol_to_line"].GetDouble()

        # Delete the previous files
        self._DeleteExistingFiles()

        # Initialize output control variables
        self.printing_times = self.settings["printing_times"].GetVector()
        self.is_printed = [False] * len(self.printing_times)

    def ExecuteBeforeSolutionLoop(self):
        """The model part related variables are initialized:
        The list of nodes and the list of variables are created.
        """
        if self.settings["use_mesh_nodes"].GetBool():
            self.nodes = []
            tolerance = self._GetTolerance()
            for node in self.model_part.Nodes:
                if self._DistanceToLine(node) < tolerance:
                    self.nodes.append(node)
        else:
            self.nodes = []
            self.elements = []
            self.area_coords = []
            configuration = KM.Configuration.Current
            locator = KM.BruteForcePointLocator(self.model_part)
            tolerance = self._GetTolerance()
            num_nodes = self.settings["number_of_nodes"].GetInt()
            for i in range(num_nodes):
                node = KM.Point(self.origin + i / (num_nodes+1) * self.direction)
                area_coords = KM.Vector()
                found_id = locator.FindElement(node, area_coords, configuration, tolerance)
                self.nodes.append(node)
                self.elements.append(self.model_part.Elements[found_id])
                self.area_coords.append(area_coords)

        self.variables, self.names = self._GenerateVariablesList(self.settings["output_variables"])
        self.nonhist_variables, self.nonhist_names = self._GenerateVariablesList(self.settings["nonhistorical_variables"], False)

    def Check(self):
        """Verify if the two specified points are found inside the domain."""
        start_point = KM.Point(self.settings["start_point"].GetVector())
        end_point = KM.Point(self.settings["end_point"].GetVector())
        locator = KM.BruteForcePointLocator(self.model_part)
        configuration = KM.Configuration.Initial
        tolerance = self._GetTolerance()
        msg = 'The {} point was not found in the domain. Please, check the geometry or the relative tolerance'
        if self.settings["use_mesh_nodes"].GetBool():
            if locator.FindNode(start_point, configuration, tolerance) < 0:
                KM.Logger.PrintWarning(msg.format('starting'))
            if locator.FindNode(end_point, configuration, tolerance) < 0:
                KM.Logger.PrintWarning(msg.format('ending'))
        else:
            area_coords = KM.Vector()
            if locator.FindElement(start_point, area_coords, configuration, tolerance) < 0:
                KM.Logger.PrintWarning(msg.format('starting'))
            if locator.FindElement(end_point, area_coords, configuration, tolerance) < 0:
                KM.Logger.PrintWarning(msg.format('ending'))

    def IsOutputStep(self):
        """This method checks if the current time step is
        near enough to the specified printing times.
        """
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        for i in range(len(self.printing_times)):
            if time >= self.printing_times[i] and not self.is_printed[i]:
                self.is_printed[i] = True
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
        if self.settings["use_mesh_nodes"].GetBool():
            for node in self.nodes:
                file.write(self._GetNodeData(node, self._DistanceToOrigin(node)))
        else:
            for node, elem, area_coords in zip(self.nodes, self.elements, self.area_coords):
                file.write(self._GetElementData(elem, area_coords, self._DistanceToOrigin(node)))
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

    def _GetNodeData(self, node, position):
        line = str(position) + " "
        for var in self.variables:
            line += str(node.GetSolutionStepValue(var)) + " "
        for var in self.nonhist_variables:
            line += str(node.GetValue(var)) + " "
        return line + "\n"

    def _GetElementData(self, elem, area_coords, position):
        line = str(position) + " "
        for var in self.variables:
            line += str(self._Interpolate(elem, area_coords, var)) + " "
        for var in self.nonhist_variables:
            line += str(self._InterpolateNonHistorical(elem, area_coords, var)) + " "
        return line + "\n"

    @staticmethod
    def _Interpolate(elem, area_coords, variable):
        nodes = elem.GetNodes()
        value = nodes[0].GetSolutionStepValue(variable) * area_coords[0]
        for n , c in zip(nodes[1:], area_coords[1:]):
            value = value + c * n.GetSolutionStepValue(variable)
        return value

    @staticmethod
    def _InterpolateNonHistorical(elem, area_coords, variable):
        nodes = elem.GetNodes()
        value = nodes[0].GetValue(variable) * area_coords[0]
        for n , c in zip(nodes[1:], area_coords[1:]):
            value = value + c * n.GetValue(variable)
        return value

    def _DeleteExistingFiles(self):
        output_path = self.settings["output_path"].GetString()
        file_name = self.settings["file_name"].GetString()
        if not output_path:
            output_path = "."
        pattern = output_path + "/" + file_name + "*"
        for filename in glob.glob(pattern):
            try:
                os.remove(filename)
            except:
                pass
