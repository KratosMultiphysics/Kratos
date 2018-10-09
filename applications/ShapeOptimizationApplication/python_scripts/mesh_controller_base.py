# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# ==============================================================================
class MeshController:
    # --------------------------------------------------------------------------
    def Initialize(self):
        pass

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, variable, design_surface):
        raise NotImplementedError("Mesh controller base class is called. Please check your implementation of the function >> UpdateMeshAccordingInputVariable << .")

# ==============================================================================
