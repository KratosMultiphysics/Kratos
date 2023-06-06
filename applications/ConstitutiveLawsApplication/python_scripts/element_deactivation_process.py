
import KratosMultiphysics as KM

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ElementDeactivationProcess(Model, settings["Parameters"])

class ElementDeactivationProcess(KM.Process):
    """This process deactivates elements according to a defined integration point variable threshold

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
        super().__init__()

        default_settings = KM.Parameters("""
        {
            "model_part_name"             : "please_specify_model_part_name",
            "variable_name"               : "DAMAGE",
            "variable_maximum_threshold"  : 0.9999,
            "average_calculation_over_ip" : true
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KM.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable_threshold = settings["variable_maximum_threshold"].GetDouble()
        self.average_variable_over_ip = settings["average_calculation_over_ip"].GetBool()

        if self.model_part.IsDistributed():
            raise RuntimeError('MPI not supported yet')

    def ExecuteFinalizeSolutionStep(self):

        if self.average_variable_over_ip: # The average > threshold
            for element in self.model_part.Elements:
                variable_elemental_data = element.CalculateOnIntegrationPoints(self.variable, self.model_part.ProcessInfo)

                average_data = 0.0
                for entry in variable_elemental_data:
                    average_data += entry
                if average_data / len(variable_elemental_data) >= self.variable_threshold:
                    element.Set(KM.ACTIVE, False)
        else: # All IP values > threshold
            for element in self.model_part.Elements:
                variable_elemental_data = element.CalculateOnIntegrationPoints(self.variable, self.model_part.ProcessInfo)
                number_failed_ip = 0

                for entry in variable_elemental_data:
                    if entry >= self.variable_threshold:
                        number_failed_ip += 1
                if number_failed_ip == len(variable_elemental_data):
                    element.Set(KM.ACTIVE, False)
