import KratosMultiphysics
from math import * # without * the test cannot be run with mathematical functions

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest


class aux_object_cpp_callback:
    def __init__(self, compiled_function):
        self.compiled_function = compiled_function

    def f(self, x, y, z, t):
        return eval(self.compiled_function)


def Factory(settings, Model):
    if not isinstance(settings,KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckScalarBaseProcess(Model, settings["Parameters"])


# All the processes python should be derived from "Process"
class CheckScalarBaseProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """This process is the base class to check analytically from a function the solution (scalar) in a certain entity belonging a certain submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"            : "This process is the base class to check analytically from a function the solution (scalar) in a certain entity belonging a certain submodelpart",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "SPECIFY_VARIABLE_NAME",
            "interval"        : [0.0, 1e30],
            "value"           : 0.0,
            "tolerance_rank"  : 3
        }
        """)

        # Copy from input
        self.settings = settings

        # Detect "End" as a tag and replace it by a large number
        if(self.settings.Has("interval")):
            if(self.settings["interval"][1].IsString() ):
                if(self.settings["interval"][1].GetString() == "End"):
                    self.settings["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" + self.settings["interval"].PrettyPrintJsonString())

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if(self.settings.Has("value")):
            if(self.settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        self.settings.ValidateAndAssignDefaults(default_settings)

        # Getting model part and interval
        self.model_part = Model[self.settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        self.interval = self.settings["interval"].GetVector()

        # Check if it is a function
        self.value_is_numeric = False
        self.is_time_function = False
        if self.settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = self.settings["value"].GetDouble()
        else:
            self.function_string = self.settings["value"].GetString()
            self.aux_function = aux_object_cpp_callback(compile(self.function_string, '', 'eval', optimize=2))

            # Must depend on time alone!
            if self.function_string.find("x") == -1 and self.function_string.find("y") == -1 and self.function_string.find("z") == -1:
                self.is_time_function = True

        # Error tolerance
        self.tolerance_rank = self.settings["tolerance_rank"].GetInt()
