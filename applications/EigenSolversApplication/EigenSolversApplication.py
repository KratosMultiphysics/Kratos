from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

IssueDeprecationWarning('EigenSolversApplication', 'please use the "LinearSolversApplication" instead')

from KratosMultiphysics.LinearSolversApplication import *
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory
