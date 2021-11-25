# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger

# ==============================================================================
class OptimizationAlgorithm:
    # --------------------------------------------------------------------------
    def __init__(self, settings, analyzer, communicator, model_part_controller):
        self.model_part_controller = model_part_controller
        self.optimization_settings = settings
        self.communicator = communicator
        self.analyzer = analyzer
        self.algorithm_settings = settings["optimization_algorithm"]
        self.design_surfaces = self.model_part_controller.GetDesignSurfaces()
        self.mappers = {}

        self.optimization_model_part = self.model_part_controller.GetOptimizationModelPart()

    # --------------------------------------------------------------------------
    def CheckApplicability( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> CheckApplicability << .")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):
        self.model_part_controller.Initialize()
        self.analyzer.InitializeBeforeOptimizationLoop()

        self._CreateMappers(self.optimization_settings["design_variables"])

        self.data_logger = data_logger.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        self.analyzer.FinalizeAfterOptimizationLoop()
        self.data_logger.FinalizeDataLogging()

    # --------------------------------------------------------------------------
    def _CreateMappers(self, mappers_settings):
        for settings in mappers_settings:
            design_surface_name = settings["model_part_name"].GetString()
            design_surface = self.design_surfaces[design_surface_name]
            settings.RemoveValue("model_part_name")

            mapper = mapper_factory.CreateMapper(design_surface, design_surface, settings["filter"])
            mapper.Initialize()
            self.mappers[design_surface_name] = mapper

    # --------------------------------------------------------------------------
    def _initializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

# ==============================================================================
