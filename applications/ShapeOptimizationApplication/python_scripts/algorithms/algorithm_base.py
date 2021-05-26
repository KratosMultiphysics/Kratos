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
    def __init__(self, settings, communicator, model_part_controller):
        pass

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
        pass
# ==============================================================================
