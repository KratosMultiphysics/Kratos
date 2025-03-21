# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import assign_scalar_variable_to_entities_process

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableToEntitiesProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorVariableToEntitiesProcess(KratosMultiphysics.Process):
    """This process assigns a given value (vector) to the entities belonging a certain submodelpart

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
            "help"                 : "This process assigns a given value (vector) to the entities belonging a certain submodelpart",
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "value"                : [0.0, 0.0, 0.0],
            "local_axes"           : {},
            "entities"             : []
        }
        """
        )

        # Detect "End" as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        if not settings.Has("value"):
            raise RuntimeError("Please specify the value to set the vector to. Example:\n" \
                               + '{\n\t"value" : [10.0, "3*t", "x+y"]\n}\n')

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if not isinstance(self.variable, KratosMultiphysics.Array1DVariable3) and not isinstance(self.variable, KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorToConditionProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.aux_processes = []


        for i_dir, var_string in enumerate(["_X", "_Y", "_Z"]):
            if not settings["value"][i_dir].IsNull():
                direction_params = KratosMultiphysics.Parameters("{}")
                direction_params.AddValue("model_part_name",settings["model_part_name"])
                direction_params.AddValue("mesh_id",settings["mesh_id"])
                direction_params.AddValue("interval",settings["interval"])
                direction_params.AddValue("value",settings["value"][i_dir])
                direction_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + var_string)
                direction_params.AddValue("local_axes",settings["local_axes"])
                direction_params.AddValue("entities",settings["entities"])
                self.aux_processes.append( assign_scalar_variable_to_entities_process.AssignScalarVariableToEntitiesProcess(Model, direction_params) )

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()
