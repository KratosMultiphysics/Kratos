from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableToConstraintProcess(Model, settings["Parameters"])

import assign_vector_variable_to_entities_process

## All the processes python should be derived from "Process"
class AssignVectorVariableToConstraintProcess(assign_vector_variable_to_entities_process.AssignVectorVariableToEntitiesProcess):
    def __init__(self, Model, settings ):
        """This process assigns a given value (vector) to the constraints belonging a certain submodelpart

        Only the member variables listed below should be accessed directly.

        Public member variables:
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        # The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                 : "This process assigns a given value (vector) to the constraints belonging a certain submodelpart",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "value"                : [10.0, 2.0, 0.0],
            "local_axes"           : {},
            "entities"             : ["constraints"]
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):
            for i in range(settings["value"].size()):
                if(settings["value"][i].IsString()):
                    raise Exception("The value can only be number for constraints:"+settings["value"][i].PrettyPrintJsonString())

        # Ensure proper entities
        if (settings["entities"].size() != 1):
            settings["entities"] = default_settings["entities"]
        else:
            if (settings["entities"][0].GetString() != "constraints"):
                settings["entities"] = default_settings["entities"]

        # Construct the base process.
        super(AssignVectorVariableToConstraintProcess, self).__init__(Model, settings)


