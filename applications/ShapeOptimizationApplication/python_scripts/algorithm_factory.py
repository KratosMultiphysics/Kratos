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
from algorithm_steepest_descent import AlgorithmSteepestDescent
from algorithm_penalized_projection import AlgorithmPenalizedProjection
import mapper_factory
import data_logger_factory

# ==============================================================================
def CreateAlgorithm( ModelPartController, Analyzer, Communicator, OptimizationSettings ):
    AlgorithmName = OptimizationSettings["optimization_algorithm"]["name"].GetString()

    Mapper = mapper_factory.CreateMapper( ModelPartController, OptimizationSettings ) 
    DataLogger = data_logger_factory.CreateDataLogger( ModelPartController, Communicator, OptimizationSettings )  

    if OptimizationSettings["optimization_algorithm"]["name"].GetString() == "steepest_descent":
        return AlgorithmSteepestDescent( ModelPartController, 
                                         Analyzer, 
                                         Communicator, 
                                         Mapper, 
                                         DataLogger, 
                                         OptimizationSettings )
    elif AlgorithmName == "penalized_projection":
        return AlgorithmPenalizedProjection( ModelPartController, 
                                             Analyzer, 
                                             Communicator, 
                                             Mapper, 
                                             DataLogger, 
                                             OptimizationSettings )
    else:
        raise NameError("The following optimization algorithm not supported by the algorithm driver (name may be misspelled): " + AlgorithmName)              

# # ==============================================================================