import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp
from pathlib import Path
import math

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SetUpPreStressedOrientedCompositeMaterials(Model, settings["Parameters"])

class ApplyPreStressingImposedStrainProcess(KM.Process):

    """This process sets the pre-stressing imposed strain according to a load factor (can depend on time) inside a defined Interval.
    This process should be called together with set_up_pre_stressed_oriented_composite_materials.py in order to work properly (it sets the initial value of imposed strain)


    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KM.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KM.Parameters(
            """
        {
            "help" : "This sets the initial conditions in terms of imposed pre-stressing strain, local axes and % participation of fiber",
            "model_part_name" : "please_specify_model_part_name",
            "interval"             : [0.0, 1e30],
            "load_factor"                : 0.001
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        if not settings.Has("load_factor"):
            raise RuntimeError("Please specify the value to set the vector to. Example:\n" + '{\n\t"load_factor" : 0.1*t*x\n}\n')

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(settings)
        self.load_factor = self.CreateFunction(settings["load_factor"])

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            print("Inside Interval!!")
            print(self.load_factor.CallFunction(0,0,0,current_time,0,0,0))



    def CreateFunction(self, value):
        function_string = ""
        if value.IsNumber():
            function_string = str(value.GetDouble())
        elif value.IsString():
            function_string = value.GetString()
        function = KM.GenericFunctionUtility(function_string)
        return function