# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.time_based_ascii_file_writer_utility as AsciiWriter
from pathlib import Path

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PrintInfoInFileProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class PrintInfoInFileProcess(KratosMultiphysics.OutputProcess):
    """This process prints a text file with the required nodal or elemental information

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class
        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        super().__init__()

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                     : "This process prints nodal/elemental information ina .txt file",
            "model_part_name"          : "please_specify_model_part_name",
            "variable_name"            : "SPECIFY_VARIABLE_NAME",
            "results_type"             : "not_provided_by_default",
            "file_name"                : "info_file.dat",
            "output_control_type"      : "step",
            "integration_point_number" : 0,
            "erase_previous_info"      : true,
            "output_interval"          : 1,
            "sum_results_from_multiple_entities" : false,
            "write_buffer_size"        : 1,
            "output_path"              : ""
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.results_type = settings["results_type"].GetString()
        if not self.results_type == "nodal_historical" and not self.results_type == "nodal_non_historical" and not self.results_type == "elemental":
            raise NameError("results_type not correct, must be nodal_historical, nodal_non_historical or elemental")
        self.is_nodal_results_type = "nodal" in self.results_type
        self.file_name = settings["file_name"].GetString()
        self.output_control_type = settings["output_control_type"].GetString()
        self.output_interval = settings["output_interval"].GetDouble()
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.instant_previous_plot = 0.0
        self.integration_point = settings["integration_point_number"].GetInt()
        self.sum_results_from_multiple_entites = settings["sum_results_from_multiple_entities"].GetBool()

        if self.model_part.IsDistributed():
            raise RuntimeError('MPI not supported yet')

        if not self.sum_results_from_multiple_entites:
            if self.is_nodal_results_type and self.model_part.NumberOfNodes() > 1:
                raise NameError("The sum_results_from_multiple_entites is false but more than one node is given...")
            if not self.is_nodal_results_type and self.model_part.NumberOfElements() > 1:
                raise NameError("The sum_results_from_multiple_entites is false but more than one element is given...")

        ascii_writer_params = KratosMultiphysics.Parameters("""{}""")
        ascii_writer_params.AddValue("file_name", settings["file_name"])
        ascii_writer_params.AddValue("output_path", settings["output_path"])
        ascii_writer_params.AddValue("write_buffer_size", settings["write_buffer_size"])
        ascii_writer_params.AddEmptyValue("file_extension")
        ascii_writer_params["file_extension"].SetString(Path(self.file_name).suffix)
        header = "# In this file we print the " + settings["results_type"].GetString() + " " + settings["variable_name"].GetString() + " in the ModelPart: " + settings["model_part_name"].GetString() + "\n\n" + "# TIME\t\t" + settings["variable_name"].GetString() + "\n"
        self.ascii_writer = AsciiWriter.TimeBasedAsciiFileWriterUtility(self.model_part, ascii_writer_params, header).file

    def PrintOutput(self):
        self.SetPreviousPlotInstant()
        if self.is_nodal_results_type:
            for node in self.model_part.Nodes:
                array_values = self.GetValueToPrint(node)
                for comp in range(len(array_values)):
                    array_values[comp] = 0.0
                break
            for node in self.model_part.Nodes:
                for comp in range(len(array_values)):
                    array_values[comp] += self.GetValueToPrint(node)[comp]
                if not self.sum_results_from_multiple_entites:
                    break
            self.PrintInFile(array_values)
        else:
            for elem in self.model_part.Elements:
                array_values = self.GetValueToPrint(elem)[self.integration_point]
                if not isinstance(array_values, (float, int)): # Can be an scalar entity...
                    for comp in range(len(array_values)):
                        array_values[comp] = 0.0
                else:
                    array_values = 0.0
                break
            for elem in self.model_part.Elements:
                if self.sum_results_from_multiple_entites:
                    for ip in range(len(self.GetValueToPrint(elem))):
                        if not isinstance(array_values, (float, int)):
                            for comp in range(len(array_values)):
                                array_values[comp] += self.GetValueToPrint(elem)[ip][comp]
                        else:
                            array_values += float(self.GetValueToPrint(elem)[ip])
                else:
                    array_values = self.GetValueToPrint(elem)[self.integration_point]
                    break
            self.PrintInFile(array_values)

    def IsOutputStep(self):
        if self.output_control_type == "step":
            return self.model_part.ProcessInfo[KratosMultiphysics.STEP] - self.instant_previous_plot >= self.output_interval
        else:
            return self.__GetTime() - self.instant_previous_plot >= self.output_interval

    def SetPreviousPlotInstant(self):
        if self.output_control_type == "step":
            self.instant_previous_plot = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        else:
            self.instant_previous_plot = self.__GetTime()

    def GetValueToPrint(self, Entity):
        if self.is_nodal_results_type:
            if self.results_type == "nodal_historical":
                return Entity.GetSolutionStepValue(self.variable)
            elif self.results_type == "nodal_non_historical":
                return Entity.GetValue(self.variable)
            else:
                raise NameError("results_type not supported")
        else:
            return Entity.CalculateOnIntegrationPoints(self.variable, self.model_part.ProcessInfo)

    def __GetTime(self):
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))

    def PrintInFile(self, values):
        self.ascii_writer.write("{0:.4e}".format(self.__GetTime()).rjust(11) + "\t")
        if not isinstance(values, (float, int)):
            for value in values:
                self.ascii_writer.write("{0:.4e}".format(value).rjust(11) + "\t")
        else:
            self.ascii_writer.write("{0:.4e}".format(values).rjust(11) + "\t")
        self.ascii_writer.write("\n")

    def ExecuteFinalize(self):
        self.ascii_writer.close()