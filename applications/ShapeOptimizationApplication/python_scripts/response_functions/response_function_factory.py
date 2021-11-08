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
from . import surface_normal_shape_change
from . import face_angle
from . import airfoil_2d_responses

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "plane_based_packaging":
        return plane_based_packaging.PlaneBasedPackaging(response_id, response_settings, model)
    elif response_type == "mesh_based_packaging":
        return mesh_based_packaging.MeshBasedPackaging(response_id, response_settings, model)
    elif response_type == "surface_normal_shape_change":
        return surface_normal_shape_change.SurfaceNormalShapeChange(response_id, response_settings, model)
    elif response_type == "face_angle":
        return face_angle.FaceAngleResponseFunction(response_id, response_settings, model)
    elif response_type == "airfoil_angle_of_attack":
        return airfoil_2d_responses.AngleOfAttackResponseFunction(response_id, response_settings, model)
    elif response_type == "airfoil_chord_length":
        return airfoil_2d_responses.ChordLengthResponseFunction(response_id, response_settings, model)
    elif response_type == "airfoil_perimeter":
        return airfoil_2d_responses.PerimeterResponseFunction(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'plane_based_packaging', 'mesh_based_packaging', 'face_angle', " +
                        "'airfoil_angle_of_attack', 'airfoil_chord_length', 'airfoil_perimeter.")
