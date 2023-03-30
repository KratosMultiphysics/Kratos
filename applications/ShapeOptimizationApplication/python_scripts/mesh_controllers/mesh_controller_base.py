# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================


# ==============================================================================
class MeshController:
    # --------------------------------------------------------------------------
    def Initialize(self):
        pass

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable(self, variable):
        raise NotImplementedError("Mesh controller base class is called. Please check your implementation of the function >> UpdateMeshAccordingInputVariable << .")

# ==============================================================================
