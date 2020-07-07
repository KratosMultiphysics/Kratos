from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

IssueDeprecationWarning('EigenSolversApplication', 'please the "LinearSolversApplication" instead')

from LinearSolversApplication import *
