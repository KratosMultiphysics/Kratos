import KratosMultiphysics as KM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.ShallowWaterApplication.postprocess.line_graph_output_process import LineGraphOutputProcess
from KratosMultiphysics.point_output_process import Interpolate


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return LineEnvelopeOutputProcess(model, settings["Parameters"])


class LineEnvelopeOutputProcess(LineGraphOutputProcess):
    """This process writes the maximum results along a line to generate a graph.
    """

    def ExecuteBeforeSolutionLoop(self):
        """Constructor of the class."""
        super().ExecuteBeforeSolutionLoop()
        self.values = [[-1e6] * len(self.variables) for _ in self.found_positions]


    def PrintOutput(self):
        """Look for the maximum value. The output file will be generated at the end of the simulation.
        This task is executed at the PrintOutput for two reasons:
        - some output values might be computed at the ExecuteBeforeOutputStep
        - this task can be skipped and executed fewer times
        """
        i = 0
        for entity, area_coords in zip(self.entities, self.area_coords):
            for v, var in enumerate(self.variables):
                value = Interpolate(var, entity, area_coords, self.historical_value)
                self.values[i][v] = max(self.values[i][v], value)
            i += 1


    def ExecuteFinalize(self):
        """The output file is created, filled and closed."""
        self.file_settings["file_name"].SetString(self.file_name)
        file = TimeBasedAsciiFileWriterUtility(self.model_part, self.file_settings, self._GetHeader()).file
        for point, var_values in zip(self.found_positions, self.values):
            file.write(self._DataToString(point, var_values))
        file.close()


    def _DataToString(self, node, values):
        data = self.print_format.format(node.X)
        data += "\t" + self.print_format.format(node.Y)
        data += "\t" + self.print_format.format(node.Z)
        for value in values:
            data += "\t" + self.print_format.format(value)
        return data + "\n"


    def _GetHeader(self):
        header = super()._GetHeader()
        header = header.replace("Results", "Envelope")
        return header
