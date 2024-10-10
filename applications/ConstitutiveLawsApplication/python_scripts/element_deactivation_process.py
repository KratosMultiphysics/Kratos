import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp

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

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.process = CLApp.ElementDeactivationProcess(self.model_part, settings)

        if self.model_part.IsDistributed():
            raise RuntimeError('MPI not supported yet')

    def ExecuteFinalizeSolutionStep(self):
        self.process.ExecuteFinalizeSolutionStep()