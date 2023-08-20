# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:
#
# ==============================================================================

import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory


class CenterlineVertexMorphingMapper():
    """
    TODO: description
    """

    def __init__(self, origin_model_part, destination_model_part, settings):

        self.settings = settings
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

        # TODO: decide if setting validation is necessary
        # in_plane_settings = self.settings["in_plane_morphing_settings"]
        # in_plane_settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultInPlaneSettings())

        self.vm_mapper = KSO.MapperCenterline(origin_model_part, destination_model_part, settings)

        self.linear_solver = None

    @classmethod
    def GetDefaultCenterlineSettings(cls):
        return KM.Parameters("""{
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : "UNKNOWN_NAME"
            },
            "background_sub_model_part_name" : "",
            "spacial_mapper_settings" : {
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }
        }""")

    def Initialize(self):

        # create linear solver
        # TODO: decide which linear solver type works best
        settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_llt" }')
        self.linear_solver = dense_linear_solver_factory.ConstructSolver(settings)

        self.vm_mapper.Initialize()

    def Update(self):

        self.vm_mapper.Update(self.linear_solver)

    def Map(self, origin_variable, destination_variable):
        self.vm_mapper.Map(origin_variable, destination_variable)

    def InverseMap(self, destination_variable, origin_variable):
        self.vm_mapper.InverseMap(destination_variable, origin_variable)
