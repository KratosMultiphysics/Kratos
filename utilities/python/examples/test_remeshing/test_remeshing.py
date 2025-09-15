import KratosMultiphysics as Kratos
from KratosMultiphysics.HDF5Application.xdmf_utils import WriteMultifileTemporalAnalysisToXdmf
from remeshing import RefineMesh
from hdf5_utilities import GetHDF5File
from hdf5_utilities import OutputModelPartToHDF5
from hdf5_utilities import OutputNodalResultsToHDF5

from math import sqrt

def write_h5(file_name, model_part):
    h5_file = GetHDF5File(file_name, "truncate")
    OutputModelPartToHDF5(model_part, h5_file)
    OutputNodalResultsToHDF5(model_part, h5_file, [Kratos.PRESSURE])
    del h5_file
    WriteMultifileTemporalAnalysisToXdmf(file_name, "/ModelData", "/ResultsData")

mmg_properties = Kratos.Parameters("""{
    "strategy"               : "hessian",
    "automatic_remesh"       : false,
    "enforce_current"        : false,
    "maximal_size"           : 100.0,
    "minimal_size"           : 1e-4,
    "hessian_strategy_parameters":{
        "metric_variable"                  : ["PRESSURE"],
        "non_historical_metric_variable"   : [false],
        "interpolation_error"              : 2e-4
    },
    "model_part_name" : "MeshingModelPart",
    "step_frequency"  : 1,
    "force_min" : true,
    "force_max" : true,
    "echo_level": 3
}""")

def SetVariableValues(model_part):
    write_h5("before_refinement.h5", model_part)

    # here set the PRESSURE variable for what ever the higher values where you want to see more refinement

    # following is an example for blocking. When Blocked, those nodes will not be changed. You can block elements and conditions as well if you like.
    # in this example, since I am not setting any PRESSURE values, it will try to remove as many nodes as possible in the un blocked regions
    # because the gradients are 0.0 everywhere. So in order to have some specialized refinement, you have to have some values in nodes so
    # that there will be some gradients.
    # if you want you can read values from h5 to get some realistic distributions as well.

    geometric_center = Kratos.Array3(0.0)
    refinement_model_part = model_part.GetSubModelPart("Slip2D.Slip2D_cylinder")
    refinement_influence_radius = 0.3
    for node in refinement_model_part.Nodes:
        geometric_center[0] += node.X
        geometric_center[1] += node.Y
        geometric_center[2] += node.Z

    geometric_center = geometric_center / len(refinement_model_part.Nodes)

    for node in model_part.Nodes:
        distance = sqrt((node.X - geometric_center[0]) ** 2 + (node.Y - geometric_center[1]) ** 2 + (node.Z - geometric_center[2]) ** 2)
        # blocking refinement if the distance is greater than refinement_influence_radius
        node.Set(Kratos.BLOCKED, distance < refinement_influence_radius)

refined_model_part = RefineMesh("cylinder", "modified_cylinder", mmg_properties, [Kratos.PRESSURE], SetVariableValues)
write_h5("after_refinement.h5", refined_model_part)