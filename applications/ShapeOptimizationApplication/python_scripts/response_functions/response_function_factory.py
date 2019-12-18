# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

from . import plane_based_packaging
from . import mesh_based_packaging

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "plane_based_packaging":
        return plane_based_packaging.PlaneBasedPackaging(response_id, response_settings, model)
    elif response_type == "mesh_based_packaging":
        return mesh_based_packaging.MeshBasedPackaging(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'plane_based_packaging', 'mesh_based_packaging'.")
