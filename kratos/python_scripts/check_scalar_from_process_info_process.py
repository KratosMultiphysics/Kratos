import KratosMultiphysics
import sys
from math import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest


class aux_object_cpp_callback:
    def __init__(self, compiled_function):
        self.compiled_function = compiled_function

    def f(self, x, y, z, t):
        return eval(self.compiled_function)


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckScalarFromProcessInfoProcess(Model, settings["Parameters"])


# All the processes python should be derived from "Process"
class CheckScalarFromProcessInfoProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """This process checks analytically from a function the solution (scalar) in the process info belonging a certain submodelpart

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
            "help"            : "This process checks analytically from a function the solution (scalar) in the process info belonging a certain submodelpart",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "SPECIFY_VARIABLE_NAME",
            "interval"        : [0.0, 1e30],
            "value"           : 0.0,
            "tolerance_rank"  : 3
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" + settings["interval"].PrettyPrintJsonString())

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        # Getting model part and interval
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = settings["interval"][0].GetDouble()
        self.interval[1] = settings["interval"][1].GetDouble()

        # Check if it is a function
        self.value_is_numeric = False
        self.is_time_function = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:
            self.function_string = settings["value"].GetString()
            if (sys.version_info > (3, 0)):
                self.aux_function = aux_object_cpp_callback(compile(self.function_string, '', 'eval', optimize=2))
            else:
                self.aux_function = aux_object_cpp_callback(compile(self.function_string, '', 'eval'))

            # Must depend on time alone!
            if self.function_string.find("x") == -1 and self.function_string.find("y") == -1 and self.function_string.find("z") == -1:
                self.is_time_function = True

            if not self.is_time_function:
                raise Exception("Must depend on time alone!: " + self.function_string)

        # Error tolerance
        self.tolerance_rank = settings["tolerance_rank"].GetInt()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        process_info = self.model_part.ProcessInfo
        current_time = process_info[KratosMultiphysics.TIME] - process_info[KratosMultiphysics.DELTA_TIME]

        if current_time >= self.interval[0] and current_time < self.interval[1]:
            if self.value_is_numeric:
                value = process_info.GetValue(self.variable)
                self.assertAlmostEqual(self.value, value, self.tolerance_rank)
            else:
                self.value = self.aux_function.f(0.0, 0.0, 0.0, current_time)
                value = process_info.GetValue(self.variable)
                self.assertAlmostEqual(self.value, value, self.tolerance_rank)
