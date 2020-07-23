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
from .packaging_response_base import PackagingResponseBase
from .. import custom_math as cm

class PlaneBasedPackaging(PackagingResponseBase):
    """
    A class that defines the response function for plane-based packaging. The plane is defined by an origin point and a normal vector.
    By default the normal of the plane indicates the feasible side of the plane (see setting 'feasible_in_normal_direction')
    """

    def __init__(self, identifier, response_settings, model):

        super().__init__(identifier, response_settings, model)

        self.plane_origin = []
        for i in range(0,3):
            self.plane_origin.append(self.response_settings["plane_origin"][i].GetDouble())

        self.plane_normal = []
        for i in range(0,3):
            self.plane_normal.append(self.response_settings["plane_normal"][i].GetDouble())

        self.plane_normal = cm.ScalarVectorProduct(1.0/cm.Norm2(self.plane_normal), self.plane_normal)

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "plane_origin"          : [0.0, 0.0, 0.0],
            "plane_normal"          : [0.0, 0.0, 1.0]
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultSettings())
        return this_defaults

    def _DistanceVectorToPlane(self, node):
        plane_to_point = [
            node.X - self.plane_origin[0],
            node.Y - self.plane_origin[1],
            node.Z - self.plane_origin[2]
        ]

        distance = cm.Dot(plane_to_point, self.plane_normal)

        return distance, self.plane_normal

    def _CalculateDistances(self):
        self.signed_distances = []
        self.directions = []
        for node in self.model_part.Nodes:
            distance, direction = self._DistanceVectorToPlane(node)
            self.signed_distances.append(distance)
            self.directions.extend(direction)
