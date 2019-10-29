from __future__ import print_function, absolute_import, division

import time as timer

import KratosMultiphysics as KM
from KratosMultiphysics import Logger

import KratosMultiphysics.ShapeOptimizationApplication as KSO

from .packaging_response_base import PackagingResponseBase
from .. import custom_math as cm

import numpy as np

class MeshPackagingResponse(PackagingResponseBase):

    def __init__(self, identifier, response_settings, model):

        super().__init__(identifier, response_settings, model)

        self.packaging_model_part = None
        self.packaging_model_part_needs_to_be_imported = False

        packaging_model_part_name = response_settings["packaging_model_part_name"].GetString()
        input_type = response_settings["packaging_model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.packaging_model_part = self.model.CreateModelPart(packaging_model_part_name, 2)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("PlanePackagingResponse: Invalid 'domain_size': {}".format(domain_size))
            self.packaging_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
            self.packaging_model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.packaging_model_part = self.model.GetModelPart(packaging_model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.packaging_model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.packaging_model_part.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)


    def Initialize(self):
        super().Initialize()

        if self.packaging_model_part_needs_to_be_imported:
            # import model part
            model_part_io = KM.ModelPartIO(self.response_settings["packaging_model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.packaging_model_part)

    def _CalculateProjectedDistances(self):

        geometry_tools = KSO.GeometryUtilities(self.model_part)

        self.signed_distances = []
        self.directions = []

        geometry_tools.ComputeProjectedDistancesToBoundingModelPart(
            self.packaging_model_part,
            self.signed_distances,
            self.directions
        )