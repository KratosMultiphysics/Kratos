from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from . import copy_one_to_many

def Create(settings):
    IssueDeprecationWarning('CopySingleToDistributed', 'please use "copy_one_to_many" instead')
    return copy_one_to_many.CopyOneToMany(settings)
