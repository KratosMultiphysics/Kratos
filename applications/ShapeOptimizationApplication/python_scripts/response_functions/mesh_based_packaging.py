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
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from .packaging_response_base import PackagingResponseBase
from ..custom_ios.wrl_io import WrlIO

class MeshBasedPackaging(PackagingResponseBase):
    """
    A class for mesh packaging response function. The mesh should contain surface conditions only.
    By default the normals of the conditions indicate the feasible side of the mesh (see setting 'feasible_in_normal_direction')
    """

    def __init__(self, identifier, response_settings, model):
        super().__init__(identifier, response_settings, model)

        self.packaging_model_part = None
        self.packaging_model_part_needs_to_be_imported = False

        packaging_model_part_name = response_settings["packaging_model_part_name"].GetString()
        self.packaging_input_type = response_settings["packaging_model_import_settings"]["input_type"].GetString()
        if self.packaging_input_type in ["mdpa", "vrml", "wrl"]:
            self.packaging_model_part = self.model.CreateModelPart(packaging_model_part_name)
            domain_size = response_settings["packaging_domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("PlaneBasedPackaging: Invalid 'domain_size': {}".format(domain_size))
            self.packaging_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
        elif self.packaging_input_type == "use_input_model_part":
            self.packaging_model_part = self.model.GetModelPart(packaging_model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.packaging_model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.packaging_model_part.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "packaging_model_part_name"       : "UNKNOWN_NAME",
            "packaging_domain_size"           : 3,
            "packaging_model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "UNKNOWN_NAME"
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultSettings())
        return this_defaults

    def Initialize(self):
        super().Initialize()

        if self.packaging_input_type == "mdpa":
            model_part_io = KM.ModelPartIO(self.response_settings["packaging_model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.packaging_model_part)
        elif self.packaging_input_type in ["vrml", "wrl"]:
            model_part_io = WrlIO(self.response_settings["packaging_model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.packaging_model_part)

    def _CalculateDistances(self):
        geometry_tools = KSO.GeometryUtilities(self.model_part)

        self.signed_distances = []
        self.directions = []

        geometry_tools.ComputeDistancesToBoundingModelPart(
            self.packaging_model_part,
            self.signed_distances,
            self.directions
        )
