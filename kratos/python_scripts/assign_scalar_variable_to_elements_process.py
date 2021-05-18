# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarVariableToElementsProcess(Model, settings["Parameters"])

from KratosMultiphysics import assign_scalar_variable_to_entities_process

## All the processes python should be derived from "Process"
class AssignScalarVariableToElementsProcess(assign_scalar_variable_to_entities_process.AssignScalarVariableToEntitiesProcess):
    """This process sets a variable a certain scalar value in a given direction, for all the entities belonging to a submodelpart. Uses assign_scalar_variable_to_elements_process for each component

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

        # The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"            : "This process assigns a given value (scalar) to the elements belonging a certain submodelpart",
            "mesh_id"         : 0,
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "SPECIFY_VARIABLE_NAME",
            "interval"        : [0.0, 1e30],
            "value"           : 0.0,
            "local_axes"      : {},
            "entities"        : ["elements"]
        }
        """
        )

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if settings.Has("value"):
            if settings["value"].IsString():
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        # Ensure proper entities
        if settings["entities"].size() != 1:
            settings["entities"] = default_settings["entities"]
        else:
            if settings["entities"][0].GetString() != "elements":
                settings["entities"] = default_settings["entities"]

        # Construct the base process.
        super(AssignScalarVariableToElementsProcess, self).__init__(Model, settings)
