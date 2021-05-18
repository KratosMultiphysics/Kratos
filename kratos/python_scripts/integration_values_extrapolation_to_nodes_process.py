# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return IntegrationValuesExtrapolationToNodesProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class IntegrationValuesExtrapolationToNodesProcess(KratosMultiphysics.Process):
    """This process extrapolates the values from integration points to the mesh nodes

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
        KratosMultiphysics.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                       : "This process extrapolates the values from integration points to the mesh nodes",
            "model_part_name"            : "",
            "echo_level"                 : 0,
            "average_variable"           : "NODAL_AREA",
            "area_average"               : true,
            "list_of_variables"          : [],
            "extrapolate_non_historical" : true
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        extrapolation_parameters = KratosMultiphysics.Parameters("""{}""")
        extrapolation_parameters.AddValue("echo_level", settings["echo_level"])
        extrapolation_parameters.AddValue("average_variable", settings["average_variable"])
        extrapolation_parameters.AddValue("area_average", settings["area_average"])
        extrapolation_parameters.AddValue("list_of_variables", settings["list_of_variables"])
        extrapolation_parameters.AddValue("extrapolate_non_historical", settings["extrapolate_non_historical"])
        self.integration_values_extrapolation_to_nodes_process = KratosMultiphysics.IntegrationValuesExtrapolationToNodesProcess(self.model_part, extrapolation_parameters)

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before the solution loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.integration_values_extrapolation_to_nodes_process.ExecuteBeforeSolutionLoop()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.integration_values_extrapolation_to_nodes_process.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        """ This method is executed at the end of the simulation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.integration_values_extrapolation_to_nodes_process.ExecuteFinalize()
