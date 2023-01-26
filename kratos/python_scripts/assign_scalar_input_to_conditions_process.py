# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarInputToConditionsProcess(Model, settings["Parameters"])

from KratosMultiphysics import assign_scalar_input_to_entities_process

## All the processes python should be derived from "Process"
class AssignScalarInputToConditionsProcess(assign_scalar_input_to_entities_process.AssignScalarInputToEntitiesProcess):
    """This process sets a variable a certain input value in a given direction, for all the conditions belonging to a submodelpart. Uses assign_scalar_input_to_conditions_process for each component

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
            "help"               : "This process assigns a given value (input) to all the conditions belonging a certain submodelpart",
            "mesh_id"            : 0,
            "model_part_name"    : "please_specify_model_part_name",
            "variable_name"      : "SPECIFY_VARIABLE_NAME",
            "interval"           : [0.0, 1e30],
            "file"               : "",
            "transfer_algorithm" : "nearest_neighbour",
            "entities"           : ["conditions"]
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        # Ensure proper entities
        if settings["entities"].size() != 1:
            settings["entities"] = default_settings["entities"]
        else:
            if settings["entities"][0].GetString() != "conditions":
                settings["entities"] = default_settings["entities"]

        # Construct the base process.
        super(AssignScalarInputToConditionsProcess, self).__init__(Model, settings)
