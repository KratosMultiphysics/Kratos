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
from KratosMultiphysics.ALEApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports

# ==============================================================================
def CreateMeshController( OptimizationModelPart, MeshSolverSettings ):
    mesh_solver_module = __import__( MeshSolverSettings["solver_type"].GetString() )
    return mesh_solver_module.CreateSolver( OptimizationModelPart, MeshSolverSettings )

# # ==============================================================================