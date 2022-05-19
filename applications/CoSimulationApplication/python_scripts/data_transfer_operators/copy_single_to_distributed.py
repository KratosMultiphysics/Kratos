from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from . import transfer_one_to_many

def Create(*args):
    IssueDeprecationWarning('CopySingleToDistributed', 'please use "transfer_one_to_many" instead')
    return transfer_one_to_many.TransferOneToMany(*args)
