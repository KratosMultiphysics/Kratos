from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateDataTransferOperator(coupling_operation_settings):
    """This function creates and returns the Data Transfer Operator used for CoSimulation"""
    return base_factory.Create(coupling_operation_settings, [], "KratosMultiphysics.CoSimulationApplication.data_transfer_operators")
