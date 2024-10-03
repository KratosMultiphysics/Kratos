# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    
    """
    This process reads the hole diameter d and the specimen width H to compute the punching stress concentration factor Kt 
    according to the equation Kt = 2 - 0284 * (1 - d/H) - 0.600 * (1 - d/H)^2 + 1.32 * (1 - d/H)^3 in the high cycle fatigue integrator.
    Reference: Pilkey WD. Peterson’s stress concentration factors. 2nd ed. John Wiley and Sons (1997).
    """

    default_settings = KM.Parameters(
        """{
            "help"                     : "This sets the inputs to compute the punching stress concentration factor Kt",
            "model_part_name"          : "please_specify_model_part_name",
            "hole_diameter"            : 0.0,
            "specimen_width"           : 1.0
        }""")
    
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    process_settings.RemoveValue("help")
    process_settings.RemoveValue("model_part_name")

    return CLA.SetPunchingStressConcentrationFactor(computing_model_part, process_settings)