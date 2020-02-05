from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignDisplacementControlProcess(Model, settings["Parameters"])

from KratosMultiphysics import assign_scalar_variable_to_entities_process as asvtep

## All the processes python should be derived from "Process"
class AssignDisplacementControlProcess(asvtep.AssignScalarVariableToEntitiesProcess):
    """This process assigns given POINT_LOAD and PRESCRIBED_DISPLACEMENT values
    (scalar) to the DisplacementControl belonging a certain submodelpart.
    Currently, it is restricted to one node and one global direction.

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

        #The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"            : "This process assigns given POINT_LOAD and PRESCRIBED_DISPLACEMENT values (scalar) to the DisplacementControl belonging a certain submodelpart",
            "model_part_name" : "please_specify_model_part_name",
            "direction"       : "x, y, or z",
            "interval"        : [0.0, 1e30],
            "point_load_value"              : 1.0,
            "prescribed_displacement_value" : "please give an expression in terms of the variable x, y, z, t"
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)
        if settings["direction"].GetString().lower() == "x":
            direction = "_X"
        elif settings["direction"].GetString().lower() == "y":
            direction = "_Y"
        elif settings["direction"].GetString().lower() == "z":
            direction = "_Z"
        else:
            raise NotImplementedError("Direction can be x, y, or z. Given is: " + settings["direction"].GetString())

        tmp = KratosMultiphysics.Parameters('{"variable_name_load": "POINT_LOAD' +direction+'", "entities": ["conditions"], "variable_name_disp" : "PRESCRIBED_DISPLACEMENT"}')

        params_load = KratosMultiphysics.Parameters("{}")
        params_load.AddValue("model_part_name", settings["model_part_name"])
        params_load.AddValue("interval", settings["interval"])
        params_load.AddValue("value", settings["point_load_value"])
        params_load.AddValue("variable_name", tmp["variable_name_load"])
        params_load.AddValue("entities", tmp["entities"])

        params_disp = KratosMultiphysics.Parameters("{}")
        params_disp.AddValue("model_part_name", settings["model_part_name"])
        params_disp.AddValue("interval", settings["interval"])
        params_disp.AddValue("value", settings["prescribed_displacement_value"])
        params_disp.AddValue("variable_name", tmp["variable_name_disp"])
        params_disp.AddValue("entities", tmp["entities"])

        # Set processes
        self.aux_processes = []
        self.aux_processes.append(asvtep.AssignScalarVariableToEntitiesProcess(Model, params_load))
        self.aux_processes.append(asvtep.AssignScalarVariableToEntitiesProcess(Model, params_disp))

        self.step_is_active = False

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()
