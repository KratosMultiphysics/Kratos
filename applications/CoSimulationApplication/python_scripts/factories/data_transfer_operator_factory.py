from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateDataTransferOperator(coupling_operation_settings, *args):
    """This function creates and returns the Data Transfer Operator used for CoSimulation"""
    return base_factory.Create(coupling_operation_settings, [*args], "KratosMultiphysics.CoSimulationApplication.data_transfer_operators")
