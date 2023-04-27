from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class MasterControl(object):

    def __init__(self) -> None:
        self.SubControlList = []
        self.actual_model_des_var
        self.Analysis = []

    def UpdateModels(self, des_var):
        for sub_control in self.SubControlList:
            sub_des_var = self.GetSubVectorPart(des_var, sub_control)
            sub_actual_model_des_var = self.GetSubVectorPart(self.actual_model_des_var)
            sub_control.UpdateModel(sub_des_var, sub_actual_model_des_var)
        self.actual_model_des_var = des_var
        for analysis in self.Analysis:
            for sub_control in self.SubControlList:
                RMP = sub_control.GetRootModel()
                if RMP.isSolved and analysis.Has(RMP):
                    analysis.isSolved = False

    def GetInitialDesignVariables(self):
        return init_des_var
    
    def FlatternVectors(self):
        vector = []
        for sub_control in self.SubControlList:
            sub_vector = sub_control.GetData()
            vector.append(sub_vector)
        return vector