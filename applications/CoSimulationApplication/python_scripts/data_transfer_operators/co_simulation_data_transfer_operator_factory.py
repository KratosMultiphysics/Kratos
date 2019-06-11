from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

def CreateDataTransferOperator(data_transfer_operator_settings):
    """
    This function creates and returns the data transfer operator.
    """
    data_transfer_operator_type = data_transfer_operator_settings["type"].GetString()
    module_full  = 'KratosMultiphysics.CoSimulationApplication.data_transfer_operators.'+data_transfer_operator_type
    module_full += '_data_transfer_operator'
    data_tranfer_operator_module = __import__(module_full, fromlist=[data_transfer_operator_type])
    return data_tranfer_operator_module.Create(data_transfer_operator_settings)
