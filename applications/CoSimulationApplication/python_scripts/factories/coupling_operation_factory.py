from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

def CreateCouplingOperation(coupling_operation_settings, solvers):
    """
    This function creates and returns the coupling operation.
    """
    coupling_operation_type = settings["type"].GetString()
    module_full  = 'KratosMultiphysics.CoSimulationApplication.coupling_operations.'+coupling_operation_type

    coupling_operation_module = __import__(module_full, fromlist=[coupling_operation_type])
    return coupling_operation_module.Create(coupling_operation_settings, solvers)
