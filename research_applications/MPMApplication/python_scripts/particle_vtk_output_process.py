import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.MPMApplication.mpm_vtk_output_process import MPMVtkOutputProcess

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    IssueDeprecationWarning("MPMApplication:","`ParticleVTKOutputProcess` is deprecated and replaced with `MPMVtkOutputProcess`")
    return MPMVtkOutputProcess(model, settings["Parameters"])
