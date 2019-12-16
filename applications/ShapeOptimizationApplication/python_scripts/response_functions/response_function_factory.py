# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

from . import plane_packaging_response
from . import mesh_packaging_response

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "plane_packaging":
        return plane_packaging_response.PlanePackagingResponse(response_id, response_settings, model)
    elif response_type == "mesh_packaging":
        return mesh_packaging_response.MeshPackagingResponse(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'plane_packaging'.")