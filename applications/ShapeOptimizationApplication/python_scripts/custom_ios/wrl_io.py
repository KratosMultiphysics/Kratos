# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

import KratosMultiphysics as KM
from .wrl_reader import read_shapes, detect_file


def _rename_to_valid_name(model_part, shape):
    name = shape.name

    ref_name = name
    i = 2
    while model_part.HasSubModelPart(name):
        name = ref_name +"_{}".format(i)
        i += 1

    if shape.name != name:
        KM.Logger.PrintWarning("ShapeOpt", "WrlIO: Name of the sub model part has been changed from '{}' to '{}' in order "\
            "to avoid name clashes!".format(shape.name, name))
        shape.name = name


class WrlIO:

    def __init__(self, file_name):
        self.file_name = detect_file(file_name)

    def ReadModelPart(self, model_part):
        KM.Logger.PrintInfo("ShapeOpt", "Start reading model part from '{}'.".format(self.file_name))

        if model_part.ProcessInfo.GetValue(KM.DOMAIN_SIZE) != 3:
            raise Exception("WrlIO: Domain size has to be 3!")

        shapes = read_shapes(self.file_name)

        nodes_shift = 0
        triangles_shift = 0
        for i, shape in enumerate(shapes):
            _rename_to_valid_name(model_part, shape)
            sub_model_part = model_part.CreateSubModelPart(shape.name)

            property_id = i+1
            new_property = model_part.CreateNewProperties(property_id)

            for i, node in enumerate(shape.nodes):
                node_id = i + nodes_shift
                new_node = model_part.CreateNewNode(node_id, *node)
                sub_model_part.AddNode(new_node, 0)

            for i, triangle in enumerate(shape.triangles):
                triangle_id = i + triangles_shift
                node_ids = [x + nodes_shift for x in triangle]

                new_condition = model_part.CreateNewCondition("SurfaceCondition3D3N", triangle_id, node_ids, new_property)
                sub_model_part.AddCondition(new_condition)

            nodes_shift += len(shape.nodes)
            triangles_shift += len(shape.triangles)
        KM.Logger.PrintInfo("ShapeOpt", "Finished reading model part.")
