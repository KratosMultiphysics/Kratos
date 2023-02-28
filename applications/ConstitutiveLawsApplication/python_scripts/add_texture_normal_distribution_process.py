import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp
import numpy as np

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return AddTextureNormalDistributionProcess(Model, settings["Parameters"])

class AddTextureNormalDistributionProcess(KM.Process):

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
            "help"                       : "This process perturbates the coordinates of the nodes according to a normal distribution",
            "model_part_name"            : "please_specify_model_part_name",
            "normal_mu"                  : "0.0",
            "normal_sigma"               : 0.1,
            "print_modified_coordinates" : true
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        if not settings.Has("load_factor"):
            raise RuntimeError("Please specify the value to set the vector to. Example:\n" + '{\n\t"load_factor" : 0.1*t*x\n}\n')

        self.model_part = Model[settings["model_part_name"].GetString()]
        # self.interval = KM.IntervalUtility(settings)
        # self.load_factor_function = self.CreateFunction(settings["load_factor"])
        self.echo = settings["echo_level"].GetInt()
        # self.old_load_factor = 1.0e30

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # current_time = self.model_part.ProcessInfo[KM.TIME]
        # step = self.model_part.ProcessInfo[KM.STEP]

        # if self.interval.IsInInterval(current_time):
        #     current_load_factor = self.load_factor_function.CallFunction(0,0,0,current_time,0,0,0)
        #     if self.echo > 0:
        #         KM.Logger.PrintInfo("ApplyPreStressingImposedStrainProcess", "Applying a load factor of " + "{0:.4e}".format(current_load_factor).rjust(11) + ".")
        #     for element in self.model_part.Elements:
        #         if element.Has(CLApp.SERIAL_PARALLEL_IMPOSED_STRAIN): # With previous imposed strain
        #             current_strain = element.GetValue(CLApp.SERIAL_PARALLEL_IMPOSED_STRAIN)
        #             if step == 1:
        #                 element.SetValue(CLApp.SERIAL_PARALLEL_IMPOSED_STRAIN, current_strain * current_load_factor)
        #             else:
        #                 element.SetValue(CLApp.SERIAL_PARALLEL_IMPOSED_STRAIN, current_strain * current_load_factor / self.old_load_factor)
        #     self.old_load_factor = current_load_factor


