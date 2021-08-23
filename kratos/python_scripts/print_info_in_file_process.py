# Importing the Kratos Library
import KratosMultiphysics

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
            "mesh_id"                  : 0,
            "model_part_name"          : "please_specify_model_part_name",
            "variable_name"            : "SPECIFY_VARIABLE_NAME",
            "variable_type"            : "nodal",
            "file_name"                : "info_file.txt",
            "output_control_type"      : "step",
            "integration_point_number" : 0,
            "erase_previous_info"      : true,
            "output_interval"          : 1
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.is_nodal_variable_type = settings["variable_type"].GetString() == "nodal"
        self.file_name = settings["file_name"].GetString()
        self.output_control_type = settings["output_control_type"].GetString()
        self.output_interval = settings["output_interval"].GetDouble()
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.instant_previous_plot = 0.0
        self.integration_point = settings["integration_point_number"].GetInt()

        open_file_aproach = "a"
        if settings["erase_previous_info"].GetBool():
            open_file_aproach = "w"

        self.plot_file = open(self.file_name, open_file_aproach)
        self.plot_file.write("# In this file we print the " + settings["variable_type"].GetString() + " " + settings["variable_name"].GetString() + " in the ModelPart: " + settings["model_part_name"].GetString() + "\n\n")
        self.plot_file.write("# TIME\t\t" + settings["variable_name"].GetString() + "\n")
        self.plot_file.close()


    def PrintOutput(self):
        self.SetPreviousPlotInstant()
        if self.is_nodal_variable_type: # In nodal values we add them
            for node in self.model_part.Nodes:
                counter = 0
                if counter == 0:
                    array_values = self.GetValueToPrint(node)
                if counter > 0:
                    for value in array_values:
                        value += self.GetValueToPrint(node)[counter]
                counter += 1
            self.plot_file = open(self.file_name, "a")
            self.plot_file.write("{0:.4e}".format(self.__GetTime()).rjust(11) + "\t")
            for value in array_values:
                self.plot_file.write("{0:.4e}".format(value).rjust(11) + "\t")
            self.plot_file.write("\n")
            self.plot_file.close()
        else: # elemental information, not adding values
            for elem in self.model_part.Elements:
                counter = 0
                if counter == 0:
                    array_values = self.GetValueToPrint(elem)
                for value in array_values:
                    if counter > 0:
                        value += self.GetValueToPrint(elem)[counter]
                    counter += 1
            self.plot_file = open(self.file_name, "a")
            self.plot_file.write("{0:.4e}".format(self.__GetTime()).rjust(11) + "\t")
            for value in array_values:
                self.plot_file.write("{0:.4e}".format(value).rjust(11) + "\t")
            self.plot_file.write("\n")
            self.plot_file.close()


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
        if self.is_nodal_variable_type:
            return Entity.GetSolutionStepValue(self.variable)
        else:
            return Entity.CalculateOnIntegrationPoints(self.variable, self.model_part.ProcessInfo)[self.integration_point]

    def __GetTime(self):
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))
