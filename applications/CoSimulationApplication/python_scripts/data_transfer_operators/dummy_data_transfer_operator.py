from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_data_transfer_operator import CoSimulationBaseDataTransferOperator

def Create(settings):
    return DummyDataTransferOperator(settings)

class DummyDataTransferOperator(CoSimulationBaseDataTransferOperator):
    pass