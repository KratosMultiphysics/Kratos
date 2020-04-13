from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos


class Formulation(ABC):
    def __init__(self, base_computing_model_part, settings):
        self.settings = settings
        self.base_computing_model_part = base_computing_model_part
        self.list_of_formulations = []
        self.move_mesh = False

    def GetName(self):
        return self.__class__.__name__

    def AddFormulation(self, formulation):
        self.list_of_formulations.append(formulation)

    def AddVariables(self):
        self.__ExecuteFormulationMethods("AddVariables")

    def AddDofs(self):
        self.__ExecuteFormulationMethods("AddDofs")

    def PrepareModelPart(self):
        self.__ExecuteFormulationMethods("PrepareModelPart")

    def Clear(self):
        self.__ExecuteFormulationMethods("Clear")

    def Check(self):
        self.__ExecuteFormulationMethods("Check")

    def Initialize(self):
        self.__ExecuteFormulationMethods("Initialize")

    def InitializeSolutionStep(self):
        self.__ExecuteFormulationMethods("InitializeSolutionStep")

    def ExecuteBeforeCouplingSolveStep(self):
        pass

    def SolveCouplingStep(self):
        for formulation in self.list_of_formulations:
            if (not formulation.IsConverged()):
                formulation.ExecuteBeforeCouplingSolveStep()
                Kratos.Logger.PrintInfo(formulation.GetName(), "Initialized formulation coupling step.")
                formulation.SolveCouplingStep()
                Kratos.Logger.PrintInfo(formulation.GetName(), "Solved formulation coupling step.")
                formulation.ExecuteAfterCouplingSolveStep()
                Kratos.Logger.PrintInfo(formulation.GetName(), "Finalized formulation coupling step.")

    def ExecuteAfterCouplingSolveStep(self):
        pass

    def FinalizeSolutionStep(self):
        self.__ExecuteFormulationMethods("FinalizeSolutionStep")

    def Finalize(self):
        self.__ExecuteFormulationMethods("Finalize")

    def GetMinimumBufferSize(self):
        min_buffer_size = 0
        for formulation in self.list_of_formulations:
            if (min_buffer_size < formulation.GetMinimumBufferSize()):
                min_buffer_size = formulation.GetMinimumBufferSize()
        return min_buffer_size

    def IsConverged(self):
        for formulation in self.list_of_formulations:
            if (not formulation.IsConverged()):
                return False
        return True

    def GetMoveMeshFlag(self):
        return self.move_mesh

    def SetCommunicator(self, communicator):
        self.__ExecuteFormulationMethods("SetCommunicator", [communicator])
        self.communicator = communicator

    def GetCommunicator(self):
        if hasattr(self, "communicator"):
            return self.communicator
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetCommunicator\" first before retrieving.")

    def IsPeriodic(self):
        if hasattr(self, "is_periodic"):
            return self.is_periodic
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetIsPeriodic\" first before checking.")

    def SetIsPeriodic(self, value):
        self.__ExecuteFormulationMethods("SetIsPeriodic", [value])
        self.is_periodic = value

    def GetBaseModelPart(self):
        return self.base_computing_model_part

    def __ExecuteFormulationMethods(self, method_name, args = []):
        for formulation in self.list_of_formulations:
            getattr(formulation, method_name)(*args)





