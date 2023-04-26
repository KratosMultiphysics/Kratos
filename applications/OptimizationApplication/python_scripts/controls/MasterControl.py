from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class MasterControl(object):

    def __init__(self) -> None:
        self.SubControlList = []
        self.actual_model_des_var

    def ComputePrimal(self, des_var):
        for sub_control in self.SubControlList:
            require_update = sub_control.CheckModelState(des_var)
            if require_update:
                sub_control.UpdateModel(des_var) # Calls VM, computes delta_x and update the sub_model part
                sub_model_part = sub_control.GetModelPart()
                sub_model_part.ComputePrimal()
                
    def ComputeAdjoint(self, des_var):
        pass