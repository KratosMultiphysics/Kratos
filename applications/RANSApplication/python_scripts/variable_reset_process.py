import KratosMultiphysics
import KratosMultiphysics.RANSApplication as KratosRANS


def Factory(settings, Model):
    if (not isinstance(Model, KratosMultiphysics.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return VariableResetProcess(Model, settings["Parameters"])


class VariableResetProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        """This process sets wall function parameters

        This process is used to activate wall function in wall conditions by setting
        RANS_IS_WALL_FUNCTION_ACTIVE to 1 in conditions. Additionally it also sets STRUCTURE and SLIP flag
        on the same nodes as well as conditions in same model part.

        Currently, all the wall conditions should be in one model part/sub model part. Having wall conditions
        in more than one model part/sub model part is not allowed.

        Args:
            Model (Kratos.Model): [description]
            settings (Kratos.Parameters): [description]

        Raises:
            Exception: [description]
        """

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name"          :"PLEASE_CHOOSE_MODEL_PART_NAME",
                "reset_variable_names_list": [],
                "reset_interval"           : 1.0,
                "historical_container"     : false
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.reset_interval = settings["reset_interval"].GetDouble()
        self.reset_variable_names_list = settings["reset_variable_names_list"].GetStringArray()
        self.historical_container = settings["historical_container"].GetBool()
        self.current_interval = 0.0

    def ExecuteInitializeSolutionStep(self):
        if (self.current_interval >= self.reset_interval):
            self.current_interval = 0.0
            for variable_name in self.reset_variable_names_list:
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                if (self.historical_container):
                    KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(variable, self.model_part.Nodes)
                    KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Resetted historical {:s} variable in {:s}.".format(variable_name, self.model_part.Name) )
                else:
                    KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(variable, self.model_part.Nodes)
                    KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Resetted non historical {:s} variable in {:s}.".format(variable_name, self.model_part.Name) )
        else:
            self.current_interval += self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

