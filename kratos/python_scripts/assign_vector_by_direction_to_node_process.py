# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorByDirectionToNodeProcess(Model, settings["Parameters"])

from KratosMultiphysics import assign_vector_by_direction_to_entity_process

## All the processes python should be derived from "Process"
class AssignVectorByDirectionToNodeProcess(assign_vector_by_direction_to_entity_process.AssignVectorByDirectionToEntityProcess):
    """This process sets a variable a certain scalar value in a given direction, for all the nodes belonging to a submodelpart. Uses assign_scalar_variable_to_nodes_process for each component

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
            "help"                 : "This process sets a variable a certain scalar value in a given direction, for all the nodes belonging to a submodelpart. Uses assign_scalar_variable_to_nodes_process for each component",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "modulus"              : 0.0,
            "direction"            : [1.0, 0.0, 0.0],
            "local_axes"           : {},
            "entities"             : ["nodes"]
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if settings.Has("modulus"):
            if settings["modulus"].IsString():
                default_settings["modulus"].SetString("0.0")
        else:
            raise RuntimeError("Please specify the modulus of the vector")

        if settings.Has("direction"):
            if settings["direction"].IsString():
                default_settings["direction"].SetString("Automatic")
        else:
            raise RuntimeError("Please specify the direction of the vector")

        # Detect "End" as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # Ensure proper entities
        if settings["entities"].size() != 1:
            settings["entities"] = default_settings["entities"]
        else:
            if settings["entities"][0].GetString() != "nodes":
                settings["entities"] = default_settings["entities"]

        # Construct the base process.
        super(AssignVectorByDirectionToNodeProcess, self).__init__(Model, settings)
