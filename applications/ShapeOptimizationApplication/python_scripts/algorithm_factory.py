# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

from algorithm_steepest_descent import AlgorithmSteepestDescent
from algorithm_penalized_projection import AlgorithmPenalizedProjection

# ==============================================================================
def CreateAlgorithm( designSurface, analyzer, mapper, communicator, optimizationSettings ):

    optimizationAlgorithm = optimizationSettings["optimization_algorithm"]["name"].GetString()

    if optimizationAlgorithm == "steepest_descent":
        return AlgorithmSteepestDescent( designSurface, analyzer, mapper, communicator, optimizationSettings )

    elif optimizationAlgorithm == "penalized_projection":
        return AlgorithmPenalizedProjection( designSurface, analyzer, mapper, communicator, optimizationSettings )  

    else:
        raise NameError("The following optimization algorithm not supported by the algorithm driver (name may be a misspelling): " + optimizationAlgorithm)              

# # ==============================================================================