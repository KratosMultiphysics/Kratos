from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorByDirectionToConditionProcess(Model, settings["Parameters"])

import assign_vector_by_direction_to_entity_process

## All the processes python should be derived from "Process"
class AssignVectorByDirectionToConditionProcess(assign_vector_by_direction_to_entity_process.AssignVectorByDirectionToEntityProcess):
    """This process sets a variable a certain scalar value in a given direction, for all the entities belonging to a submodelpart. Uses assign_scalar_variable_to_conditions_process for each component

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, Model, settings ):
        # The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                 : "This process sets a variable a certain scalar value in a given direction, for all the conditions belonging to a submodelpart. Uses assign_scalar_variable_to_conditions_process for each component",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "modulus"              : 1.0,
            "direction"            : [1.0, 0.0, 0.0],
            "local_axes"           : {},
            "entities"             : ["conditions"]
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("modulus")):
            if(settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")

        if(settings.Has("direction")):
            if(settings["direction"].IsString()):
                default_settings["direction"].SetString("Automatic")

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # Ensure proper entities
        if (settings["entities"].size() != 1):
            settings["entities"] = default_settings["entities"]
        else:
            if (settings["entities"][0].GetString() != "conditions"):
                settings["entities"] = default_settings["entities"]

        # Construct the base process.
        super(AssignVectorByDirectionToConditionProcess, self).__init__(Model, settings)
