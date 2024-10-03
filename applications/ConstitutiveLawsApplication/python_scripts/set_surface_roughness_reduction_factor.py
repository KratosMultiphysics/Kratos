# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    
    """
    This process reads the surface roughness r_z and material constants a_r and sigma_uts_min to compute the roughness reduction factor K_r 
    according to the equation K_r = 1 - a_r * log10(r_z) * log10(2 * sigma_uts / sigma_uts_min) in the high cycle fatigue integrator.
    For steel: a_r = 0.22 and sigma_uts_min = 400 MPa
    For aluminum: a_r = 0.22 and sigma_uts_min = 133 MPa
    Reference: E. Haibach, Fkm-guideline: Analytical strength assessment of components in mechanical engineering, 5 revised edition, english version (2003).
    """

    default_settings = KM.Parameters(
        """{
            "help"                     : "This sets the inputs to compute the fatigue surface roughness reduction factor Kr",
            "model_part_name"          : "please_specify_model_part_name",
            "surface_roughness"        : 1.0,
            "mat_param_ar"             : 0.0,
            "mat_param_sigmautsmin"    : 1.0
        }""")
    
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    process_settings.RemoveValue("help")
    process_settings.RemoveValue("model_part_name")

    return CLA.SetSurfaceRoughnessReductionFactor(computing_model_part, process_settings)