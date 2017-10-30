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
def CreateAlgorithm( InputModelPart, Analyzer, OptimizationSettings ):
    AlgorithmName = OptimizationSettings["optimization_algorithm"]["name"].GetString()

    DesignSurface = GetDesignSurfaceFromInputModelPart( InputModelPart, OptimizationSettings )
    DampingRegions = GetDampingRegionsFromInputModelPart( InputModelPart, OptimizationSettings )

    Mapper = mapper_factory.CreateMapper( DesignSurface, OptimizationSettings ) 
    Communicator = communicator_factory.CreateCommunicator( OptimizationSettings )
    DataLogger = optimization_data_logger_factory.CreateDataLogger( InputModelPart, Communicator, OptimizationSettings )

    if AlgorithmName == "steepest_descent":
        return AlgorithmSteepestDescent( DesignSurface, DampingRegions, Analyzer, Mapper, Communicator, DataLogger, OptimizationSettings )
    elif AlgorithmName == "penalized_projection":
        return AlgorithmPenalizedProjection( DesignSurface, DampingRegions, Analyzer, Mapper, Communicator, DataLogger, OptimizationSettings )  
    else:
        raise NameError("The following optimization algorithm not supported by the algorithm driver (name may be misspelled): " + AlgorithmName)              

# --------------------------------------------------------------------------
def GetDesignSurfaceFromInputModelPart( InputModelPart, OptimizationSettings ):
    nameOfDesingSurface = OptimizationSettings["design_variables"]["design_submodel_part_name"].GetString()
    if InputModelPart.HasSubModelPart( nameOfDesingSurface ):
        optimizationModel = InputModelPart.GetSubModelPart( nameOfDesingSurface )
        print("> The following design surface was defined:\n\n",optimizationModel)
        return optimizationModel
    else:
        raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesingSurface)         

# --------------------------------------------------------------------------
def GetDampingRegionsFromInputModelPart( InputModelPart, OptimizationSettings ):
    DampingRegions = {}
    if(OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()):
        print("> The following damping regions are defined: \n")
        for regionNumber in range(OptimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
            regionName = OptimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
            if InputModelPart.HasSubModelPart(regionName):
                print(regionName)
                DampingRegions[regionName] = InputModelPart.GetSubModelPart(regionName)
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)    
        print("")    
    return DampingRegions   

# # ==============================================================================