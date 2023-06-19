import KratosMultiphysics as KM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.ShallowWaterApplication.postprocess.line_graph_output_process import LineGraphOutputProcess
from KratosMultiphysics.point_output_process import Interpolate


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object")
    return LineEnvelopeOutputProcess(model, settings["Parameters"])


class LineEnvelopeOutputProcess(LineGraphOutputProcess):
    """This process writes the maximum results along a line to generate a graph."""

    def ExecuteBeforeSolutionLoop(self):
        """Initialize the list of maximum values."""
        super().ExecuteBeforeSolutionLoop()
        num_of_vaviables = len(self.variables) + len(self.nonhistorical_variables)
        self.values = [[-1e6] * num_of_vaviables for _ in self.found_positions]


    def ExecuteFinalizeSolutionStep(self):
        """Look for the maximum value."""
        i = 0
        for entity, area_coords in zip(self.entities, self.area_coords):
            for v, var in enumerate(self.variables):
                value = Interpolate(var, entity, area_coords, historical_value=True)
                self.values[i][v] = max(self.values[i][v], value)
            for v, var in enumerate(self.nonhistorical_variables):
                value = Interpolate(var, entity, area_coords, historical_value=False)
                self.values[i][v] = max(self.values[i][v], value)
            i += 1


    def PrintOutput(self):
        """The output file is created, filled and closed.

        The previous output files are overwitten. If the simulation
        does not reach the end, an envelope will be kept.
        """
        self.file_settings["file_name"].SetString(self.file_name)
        file = TimeBasedAsciiFileWriterUtility(self.model_part, self.file_settings, self._GetHeader()).file
        for point, var_values in zip(self.found_positions, self.values):
            file.write(self._DataToString(point, var_values))
        file.close()


    def _DataToString(self, node, values):
        data = self.print_format.format(node.X)
        data += " " + self.print_format.format(node.Y)
        data += " " + self.print_format.format(node.Z)
        for value in values:
            data += " " + self.print_format.format(value)
        return data + "\n"


    def _GetHeader(self):
        header = super()._GetHeader()
        header = header.replace("Results", "Envelope")
        return header
