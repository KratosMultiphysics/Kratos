from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateDataTransferOperator(coupling_operation_settings):
    """This function creates and returns the Data Transfer Operator used for CoSimulation"""
    return base_factory.Create(coupling_operation_settings, [], "KratosMultiphysics.CoSimulationApplication.data_transfer_operators")
