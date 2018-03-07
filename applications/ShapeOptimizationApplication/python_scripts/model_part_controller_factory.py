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

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import time as timer

# ==============================================================================
def CreateController( OptimizationModelPart, OptimizationSettings ):
    return ModelPartController( OptimizationModelPart, OptimizationSettings )

# ==============================================================================
class ModelPartController:
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.OptimizationSettings = OptimizationSettings

        self.DesignSurface = None
        self.DampingRegions = None

        MeshMotionSettings = self.OptimizationSettings["design_variables"]["mesh_motion"]
        if MeshMotionSettings["apply_ale_mesh_solver"].GetBool():
            from mesh_controller_ale_solver import MeshControllerUsingALESolver
            self.MeshController = MeshControllerUsingALESolver( self.OptimizationModelPart, MeshMotionSettings)
        else:
            from mesh_controller_basic_updating import MeshControllerBasicUpdating
            self.MeshController = MeshControllerBasicUpdating( self.OptimizationModelPart )

    # --------------------------------------------------------------------------
    def IsOptimizationModelPartAlreadyImported( self ):
        if self.OptimizationModelPart.NumberOfNodes()>0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def ImportOptimizationModelPart( self ):
        model_part_io = ModelPartIO( self.OptimizationSettings["design_variables"]["optimization_model_part_name"].GetString() )
        model_part_io.ReadModelPart( self.OptimizationModelPart )
        self.OptimizationModelPart.SetBufferSize( 1 )
        self.OptimizationModelPart.ProcessInfo.SetValue( DOMAIN_SIZE, self.OptimizationSettings["design_variables"]["domain_size"].GetInt() )

    # --------------------------------------------------------------------------
    def InitializeMeshController( self ):
        self.MeshController.Initialize()

    # --------------------------------------------------------------------------
    def CloneTimeStep( self, new_step ):
        self.OptimizationModelPart.CloneTimeStep( new_step )

    # --------------------------------------------------------------------------
    def UpdateMeshAccordingInputVariable( self, InputVariable ):
        self.MeshController.UpdateMeshAccordingInputVariable( InputVariable )

    # --------------------------------------------------------------------------
    def SetMeshToReferenceMesh( self ):
        MeshControllerUtilities( self.OptimizationModelPart ).SetMeshToReferenceMesh()

    # --------------------------------------------------------------------------
    def SetReferenceMeshToMesh( self ):
        MeshControllerUtilities( self.OptimizationModelPart ).SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def SetDeformationVariablesToZero( self ):
        MeshControllerUtilities( self.OptimizationModelPart ).SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def GetOptimizationModelPart( self ):
        return self.OptimizationModelPart

    # --------------------------------------------------------------------------
    def GetDesignSurface( self ):
        if self.DesignSurface is None:
            self.__IdentifyDesignSurface()
        return self.DesignSurface

    # --------------------------------------------------------------------------
    def GetDampingRegions( self ):
        if self.DampingRegions is None:
            self.__IdentifyDampingRegions()
        return self.DampingRegions

    # --------------------------------------------------------------------------
    def __IdentifyDesignSurface( self ):
        nameOfDesingSurface = self.OptimizationSettings["design_variables"]["design_surface_sub_model_part_name"].GetString()
        if self.OptimizationModelPart.HasSubModelPart( nameOfDesingSurface ):
            self.DesignSurface = self.OptimizationModelPart.GetSubModelPart( nameOfDesingSurface )
            print("> The following design surface was defined:\n\n",self.DesignSurface)
        else:
            raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesingSurface)

    # --------------------------------------------------------------------------
    def __IdentifyDampingRegions( self ):
        print("> The following damping regions are defined: \n")
        self.DampingRegions = {}
        if self.OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool():
            if self.OptimizationSettings["design_variables"]["damping"].Has("damping_regions"):
                for regionNumber in range(self.OptimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
                    regionName = self.OptimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
                    if self.OptimizationModelPart.HasSubModelPart(regionName):
                        print(regionName)
                        self.DampingRegions[regionName] = self.OptimizationModelPart.GetSubModelPart(regionName)
                    else:
                        raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)
            else:
                raise ValueError("Definition of damping regions required but not availabe!")
        print("")

# ==============================================================================