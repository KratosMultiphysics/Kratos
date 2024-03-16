import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.MPMApplication.mpm_gid_output_process import MPMGiDOutputProcess

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a model object")
    wrng_msg = "Process `particle_gid_output_process` is replaced by `mpm_gid_output_process`."
    IssueDeprecationWarning("MPMApplication:",wrng_msg)
    model_part = model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]
    return ParticleGiDOutputProcess(model_part, output_name, postprocess_parameters)
