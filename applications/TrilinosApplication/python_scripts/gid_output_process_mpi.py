
import KratosMultiphysics as KM

from KratosMultiphysics.mpi.distributed_gid_output_process import DistributedGiDOutputProcess
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]

    IssueDeprecationWarning("GiDOutputProcessMPI",
        "Attempting to create deprecated process \"GiDOutputProcessMPI\",",
        "please use \"DistributedGiDOutputProcess\" in module KratosMultiphysics.mpi instead.")
    return DistributedGiDOutputProcess(model_part, output_name, postprocess_parameters)
