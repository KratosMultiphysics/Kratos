# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import assign_scalar_variable_process

from math import *

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorVariableProcess(KratosMultiphysics.Process):
    """This process assigns a given value (vector) to the nodes belonging a certain submodelpart

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

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                 : "This process assigns a given value (vector) to the nodes belonging a certain submodelpart",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "value"                : [10.0, "3*t", "x+y"],
            "constrained"          : [true,true,true],
            "local_axes"           : {}
        }
        """
        )
        #example of admissible values for "value" : [10.0, "3*t", "x+y"]

        ## Trick to ensure that if someone sets constrained as a single bool, it is transformed to a vector
        if settings.Has("constrained"):
            if settings["constrained"].IsBool():
                is_fixed = settings["constrained"].GetBool()
                #print("is_fixed = ",is_fixed)
                settings["constrained"] = default_settings["constrained"]
                for i in range(3):
                    settings["constrained"][i].SetBool(is_fixed)

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if not isinstance(self.variable, KratosMultiphysics.Array1DVariable3) and not isinstance(self.variable, KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.aux_processes = []


        # Loop over components X, Y and Z
        for indice,variable in enumerate(["_X", "_Y", "_Z"]):
            if not settings["value"][indice].IsNull():
                i_params = KratosMultiphysics.Parameters("{}")
                i_params.AddValue("model_part_name",settings["model_part_name"])
                i_params.AddValue("mesh_id",settings["mesh_id"])
                i_params.AddEmptyValue("constrained").SetBool(settings["constrained"][indice].GetBool())
                i_params.AddValue("interval",settings["interval"])
                i_params.AddValue("value",settings["value"][indice])
                i_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + variable)
                i_params.AddValue("local_axes",settings["local_axes"])
                self.aux_processes.append( assign_scalar_variable_process.AssignScalarVariableProcess(Model, i_params) )

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed in before initialize the solution step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.ExecuteInitializeSolutionStep()

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
