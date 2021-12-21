from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from . import copy_single_to_multiple

def Create(settings):
    IssueDeprecationWarning('CopySingleToDistributed', 'please use "copy_single_to_multiple" instead')
    return copy_single_to_multiple.CopySingleToMultiple(settings)
