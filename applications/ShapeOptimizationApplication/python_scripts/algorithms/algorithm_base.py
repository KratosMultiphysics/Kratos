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
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory

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
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> InitializeOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> FinalizeOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def _CreateMappers(self, mappers_settings):
        for settings in mappers_settings:
            design_surface_name = settings["design_surface_name"].GetString()
            design_surface = self.design_surfaces[design_surface_name]
            settings.RemoveValue("design_surface_name")

            mapper = mapper_factory.CreateMapper(design_surface, design_surface, settings)
            mapper.Initialize()
            self.mappers[design_surface_name] = mapper


# ==============================================================================
