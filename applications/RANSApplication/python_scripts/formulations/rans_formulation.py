import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

class RansFormulation:
    def __init__(self, base_computing_model_part, settings):
        self.settings = settings
        self.base_computing_model_part = base_computing_model_part
        self.list_of_formulations = []
        self.list_of_processes = []
        self.move_mesh = False

    def GetName(self):
        return self.__class__.__name__

    def AddRansFormulation(self, formulation):
        if (isinstance(formulation, RansFormulation)):
            self.list_of_formulations.append(formulation)
        else:
            msg = str(formulation).rstrip() + " is not a RansFormulation. Please use only RansFormulation objects."
            raise Exception(msg)

    def AddProcess(self, process):
        if (isinstance(process, KratosRANS.RansFormulationProcess)):
            self.list_of_processes.append(process)
        else:
            msg = str(process).rstrip() + " is not a RansFormulationProcess. Please use only RansFormulationProcess objects."
            raise Exception(msg)

    def AddVariables(self):
        self.__ExecuteRansFormulationMethods("AddVariables")

    def AddDofs(self):
        self.__ExecuteRansFormulationMethods("AddDofs")

    def PrepareModelPart(self):
        self.__ExecuteRansFormulationMethods("PrepareModelPart")

    def Clear(self):
        self.__ExecuteRansFormulationMethods("Clear")

    def Check(self):
        self.__ExecuteProcessMethods("Check")
        self.__ExecuteRansFormulationMethods("Check")

    def Initialize(self):
        self.__ExecuteProcessMethods("ExecuteInitialize")
        self.__ExecuteRansFormulationMethods("Initialize")

    def InitializeSolutionStep(self):
        self.__ExecuteProcessMethods("ExecuteInitializeSolutionStep")
        self.__ExecuteRansFormulationMethods("InitializeSolutionStep")

    def SolveCouplingStep(self):
        max_iterations = self.GetMaxCouplingIterations()
        for iteration in range(max_iterations):
            self.ExecuteBeforeCouplingSolveStep()
            for formulation in self.list_of_formulations:
                if (not formulation.SolveCouplingStep()):
                    return False
            self.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved coupling itr. " + str(iteration + 1) + "/" + str(max_iterations) + ".")

        return True

    def ExecuteBeforeCouplingSolveStep(self):
        self.__ExecuteProcessMethods("ExecuteBeforeCouplingSolveStep")

    def ExecuteAfterCouplingSolveStep(self):
        self.__ExecuteProcessMethods("ExecuteAfterCouplingSolveStep")

    def FinalizeSolutionStep(self):
        self.__ExecuteRansFormulationMethods("FinalizeSolutionStep")
        self.__ExecuteProcessMethods("ExecuteFinalizeSolutionStep")

    def Finalize(self):
        self.__ExecuteRansFormulationMethods("Finalize")
        self.__ExecuteProcessMethods("ExecuteFinalize")

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

        if (self.GetStrategy() is not None):
            is_converged = self.GetStrategy().IsConverged()
            if (is_converged):
                Kratos.Logger.PrintInfo(self.GetName(), " *** CONVERGENCE ACHIEVED ***")
            return is_converged

        return True

    def GetMoveMeshFlag(self):
        return self.move_mesh

    def SetCommunicator(self, communicator):
        self.__ExecuteRansFormulationMethods("SetCommunicator", [communicator])
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
        self.__ExecuteRansFormulationMethods("SetIsPeriodic", [value])
        self.is_periodic = value

    def SetTimeSchemeSettings(self, settings):
        self.__ExecuteRansFormulationMethods("SetTimeSchemeSettings", [settings])
        self.time_scheme_settings = settings

    def SetWallFunctionSettings(self, settings):
        self.__ExecuteRansFormulationMethods("SetWallFunctionSettings", [settings])
        self.wall_function_settings = settings

    def GetTimeSchemeSettings(self):
        if (hasattr(self, "time_scheme_settings")):
            return self.time_scheme_settings
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetTimeSchemeSettings\" first before calling \"GetTimeSchemeSettings\".")

    def GetBaseModelPart(self):
        return self.base_computing_model_part

    def SetMaxCouplingIterations(self, max_iterations):
        self.max_coupling_iterations = max_iterations

    def GetMaxCouplingIterations(self):
        if (hasattr(self, "max_coupling_iterations")):
            return self.max_coupling_iterations
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetMaxCouplingIterations\" first before calling \"GetMaxCouplingIterations\".")

    def SetConstants(self, settings):
        self.__ExecuteRansFormulationMethods("SetTimeSchemeSettings", [settings])

    def GetRansFormulationsList(self):
        return self.list_of_formulations

    def GetProcessList(self):
        return self.list_of_processes

    def GetModelPart(self):
        return None

    def GetStrategy(self):
        return None

    def GetInfo(self):
        info = "\n" + self.GetName()
        if (self.GetModelPart() is not None):
            info += "\n   Model part    : " + str(self.GetModelPart().Name)

        if (str(self.GetMaxCouplingIterations()) != "N/A"):
            info += "\n   Max iterations: " + str(self.GetMaxCouplingIterations())

        if (len(self.GetProcessList()) > 0):
            info += "\n   Process list:"
            for process in self.GetProcessList():
                info += "\n      " + str(process).strip()

        for formulation in self.list_of_formulations:
            info += str(formulation.GetInfo()).replace("\n", "\n   ")
        return info

    def __ExecuteRansFormulationMethods(self, method_name, args = []):
        for formulation in self.list_of_formulations:
            getattr(formulation, method_name)(*args)

    def __ExecuteProcessMethods(self, method_name):
        for process in self.list_of_processes:
            getattr(process, method_name)()




