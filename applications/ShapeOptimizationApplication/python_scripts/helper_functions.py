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

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def GetDesignSurfaceFromOptimizationModelPart( OptimizationModelPart, OptimizationSettings ):
    nameOfDesingSurface = OptimizationSettings["design_variables"]["design_surface_sub_model_part_name"].GetString()
    if OptimizationModelPart.HasSubModelPart( nameOfDesingSurface ):
        DesignSurface = OptimizationModelPart.GetSubModelPart( nameOfDesingSurface )
        print("> The following design surface was defined:\n\n",DesignSurface)
        return DesignSurface
    else:
        raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesingSurface)         

# --------------------------------------------------------------------------
def GetDampingRegionsFromOptimizationModelPart( OptimizationModelPart, OptimizationSettings ):
    print("> The following damping regions are defined: \n")
    DampingRegions = {}
    print(OptimizationSettings)
    if OptimizationSettings["design_variables"]["damping"].Has("damping_regions"):
        for regionNumber in range(OptimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
            regionName = OptimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
            if OptimizationModelPart.HasSubModelPart(regionName):
                print(regionName)
                DampingRegions[regionName] = OptimizationModelPart.GetSubModelPart(regionName)
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)  
    else:
                raise ValueError("Definition of damping regions required but not availabe!")  
    print("")    
    return DampingRegions   

# # ==============================================================================