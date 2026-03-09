import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.MPMApplication.mpm_from_json_check_result_process import MPMFromJsonCheckResultProcess

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a model object")
    wrng_msg = "Process `particle_from_json_check_result_process` is replaced by `mpm_from_json_check_result_process`."
    IssueDeprecationWarning("MPMApplication:",wrng_msg)
    return MPMFromJsonCheckResultProcess(model, settings["Parameters"])
