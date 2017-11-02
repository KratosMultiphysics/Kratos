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
import mapper_factory
import communicator_factory
import optimization_data_logger_factory as optimization_data_logger_factory

from algorithm_steepest_descent import AlgorithmSteepestDescent
from algorithm_penalized_projection import AlgorithmPenalizedProjection

# ==============================================================================
def CreateAlgorithm( OptimizationModelPart, Analyzer, OptimizationSettings ):
    AlgorithmName = OptimizationSettings["optimization_algorithm"]["name"].GetString()

    DesignSurface = GetDesignSurfaceFromOptimizationModelPart( OptimizationModelPart, OptimizationSettings )
    DampingRegions = GetDampingRegionsFromOptimizationModelPart( OptimizationModelPart, OptimizationSettings )

    Mapper = mapper_factory.CreateMapper( DesignSurface, OptimizationSettings ) 
    Communicator = communicator_factory.CreateCommunicator( OptimizationSettings )
    DataLogger = optimization_data_logger_factory.CreateDataLogger( OptimizationModelPart, DesignSurface, Communicator, OptimizationSettings )

    if AlgorithmName == "steepest_descent":
        return AlgorithmSteepestDescent( DesignSurface, DampingRegions, Analyzer, Mapper, Communicator, DataLogger, OptimizationSettings )
    elif AlgorithmName == "penalized_projection":
        return AlgorithmPenalizedProjection( DesignSurface, DampingRegions, Analyzer, Mapper, Communicator, DataLogger, OptimizationSettings )  
    else:
        raise NameError("The following optimization algorithm not supported by the algorithm driver (name may be misspelled): " + AlgorithmName)              

# --------------------------------------------------------------------------
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
    DampingRegions = {}
    if(OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()):
        print("> The following damping regions are defined: \n")
        for regionNumber in range(OptimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
            regionName = OptimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
            if OptimizationModelPart.HasSubModelPart(regionName):
                print(regionName)
                DampingRegions[regionName] = OptimizationModelPart.GetSubModelPart(regionName)
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)    
        print("")    
    return DampingRegions   

# # ==============================================================================