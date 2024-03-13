import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.MPMApplication.assign_gravity_to_material_point_process import AssignGravityToMaterialPointProcess

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a model object")
    wrng_msg  = "Process `assign_gravity_to_particle_process` is replaced "
    wrng_msg += "by `assign_gravity_to_material_point_process`."
    IssueDeprecationWarning("MPMApplication:",wrng_msg)
    return AssignGravityToMaterialPointProcess(model, settings["Parameters"])
