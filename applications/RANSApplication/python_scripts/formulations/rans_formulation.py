from abc import abstractmethod, ABC

import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.process_factory import KratosProcessFactory

class RansFormulation(ABC):
    def __init__(self, base_computing_model_part, settings):
        """RansFormulation base class

        This class is the base class for formulations used in RANSApplication. A single leaf formulation
        is responsible for solving for one variable only. RansFormulations can be added to another RansFormulation
        creating coupled hierarchical formulation. Then this base RansFormulations solves them recursively according
        to the addition order.

        Args:
            base_computing_model_part (Kratos.ModelPart): Base model part, which is copied to create solvers for given formulation
            settings (Kratos.Parameters): Settings to be used in this formulation
        """
        self.__settings = settings
        self.__base_computing_model_part = base_computing_model_part
        self.__list_of_formulations = []
        self.__list_of_processes = []
        self.__move_mesh = False

    @abstractmethod
    def GetDefaultParameters(self):
        """Returns default parameters used in this formulation

        Returns:
            Kratos.Parameters: Parameters of this formulation
        """
        pass

    def GetParameters(self):
        """Returns parameters used in this formulation

        Returns:
            Kratos.Parameters: Parameters of this formulation
        """
        return self.__settings

    def BackwardCompatibilityHelper(self, settings, deprecated_settings_dict):
        """Recursively calls BackwardCompatibilityHelper methods of existing formulations in this formulaton
        """
        self.__ExecuteRansFormulationMethods("BackwardCompatibilityHelper", [settings, deprecated_settings_dict])

    def GetDomainSize(self):
        """Returns domain size
        """
        return self.__base_computing_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

    def AddRansFormulation(self, formulation):
        """Adds another RansFormulation to current formulation creating a list of formulations.

        Args:
            formulation (RansFormulation): Formulation to be added
        """
        if (isinstance(formulation, RansFormulation)):
            self.__list_of_formulations.append(formulation)
        else:
            msg = str(formulation).rstrip() + " is not a RansFormulation. Please use only RansFormulation objects."
            raise Exception(msg)

    def AddProcess(self, process):
        """Adds a RansFormulationProcess to current RansFormulation

        Args:
            process (Kratos.RANSApplication.RansFormulationProcess): RansFormulationProcess to be added to current formulation
        """
        if (isinstance(process, KratosRANS.RansFormulationProcess)):
            self.__list_of_processes.append(process)
        else:
            msg = str(process).rstrip() + " is not a RansFormulationProcess. Please use only RansFormulationProcess objects."
            raise Exception(msg)

    def AddProcessesList(self, kratos_parameters_processes_list):
        factory = KratosProcessFactory(self.GetBaseModelPart().GetModel())
        self.auxiliar_process_list = factory.ConstructListOfProcesses(kratos_parameters_processes_list)
        for process in self.auxiliar_process_list:
            self.AddProcess(process)

    def AddVariables(self):
        """Recursively calls AddVariables methods of existing formulations in this formulaton
        """
        self.__ExecuteRansFormulationMethods("AddVariables")

    def AddDofs(self):
        """Recursively calls AddDofs methods of existing formulations in this formulaton
        """
        self.__ExecuteRansFormulationMethods("AddDofs")

    def PrepareModelPart(self):
        """Recursively calls PrepareModelPart methods of existing formulations in this formulaton
        """
        self.__ExecuteRansFormulationMethods("PrepareModelPart")

    def Clear(self):
        """Recursively calls Clear methods of existing formulations in this formulaton and clears strategy
        """
        if (self.GetStrategy() is not None):
            self.GetStrategy().Clear()

        self.__ExecuteRansFormulationMethods("Clear")

    def Check(self):
        """Recursively calls Check methods of existing formulations, processes and strategy in this formulaton
        """
        self.__ExecuteProcessMethods("Check")
        self.__ExecuteRansFormulationMethods("Check")

        if (self.GetStrategy() is not None):
            self.GetStrategy().Check()

    def Initialize(self):
        """Recursively calls Initialize methods of existing formulations, processes and strategy in this formulaton
        """
        self.__ExecuteProcessMethods("ExecuteInitialize")
        self.__ExecuteRansFormulationMethods("Initialize")

        if (self.GetStrategy() is not None):
            self.GetStrategy().Initialize()

    def InitializeSolutionStep(self):
        """Recursively calls InitializeSolutionStep methods of existing formulations, processes and strategy in this formulaton
        """
        self.__ExecuteProcessMethods("ExecuteInitializeSolutionStep")
        self.__ExecuteRansFormulationMethods("InitializeSolutionStep")

        if (self.GetStrategy() is not None):
            self.GetStrategy().InitializeSolutionStep()

    def Predict(self):
        for formulation in self.__list_of_formulations:
            formulation.Predict()

        if (self.GetStrategy() is not None):
            self.GetStrategy().Predict()

    def SolveCouplingStep(self):
        """Solves current formulation

        This method recursively solves each formulation in the list of formulations.

        Returns:
            bool: True if solve is successfull, False if not
        """
        max_iterations = self.GetMaxCouplingIterations()
        for iteration in range(max_iterations):
            self.ExecuteBeforeCouplingSolveStep()
            for formulation in self.__list_of_formulations:
                if (not formulation.SolveCouplingStep()):
                    return False
            self.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved coupling itr. " + str(iteration + 1) + "/" + str(max_iterations) + ".")
            if (self.IsConverged()):
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Coupling iter. {:d}/{:d} - *** CONVERGENCE ACHIEVED ***".format(iteration + 1, max_iterations))
                return True

        return True

    def ExecuteBeforeCouplingSolveStep(self):
        """Recursively calls ExecuteBeforeCouplingSolveStep methods of existing formulations in this formulaton
        """
        self.__ExecuteProcessMethods("ExecuteBeforeCouplingSolveStep")

    def ExecuteAfterCouplingSolveStep(self):
        """Recursively calls ExecuteAfterCouplingSolveStep methods of existing formulations in this formulaton
        """
        self.__ExecuteProcessMethods("ExecuteAfterCouplingSolveStep")

    def FinalizeSolutionStep(self):
        """Recursively calls FinalizeSolutionStep methods of existing formulations, processes and strategy in this formulaton
        """
        if (self.GetStrategy() is not None):
            self.GetStrategy().FinalizeSolutionStep()

        self.__ExecuteRansFormulationMethods("FinalizeSolutionStep")
        self.__ExecuteProcessMethods("ExecuteFinalizeSolutionStep")

    def Finalize(self):
        """Recursively calls Finalize methods of existing formulations and processes in this formulaton
        """
        self.__ExecuteRansFormulationMethods("Finalize")
        self.__ExecuteProcessMethods("ExecuteFinalize")

        if (self.GetStrategy() is not None):
            self.GetStrategy().Clear()

    def GetMinimumBufferSize(self):
        """Recursively calculate minimum buffer size required by all formulations

        Returns:
            int: Minimum buffer size
        """
        min_buffer_size = 0
        for formulation in self.__list_of_formulations:
            if (min_buffer_size < formulation.GetMinimumBufferSize()):
                min_buffer_size = formulation.GetMinimumBufferSize()
        return min_buffer_size

    def IsBufferInitialized(self):
        """Check whether enough buffer is initialized to solve current formulation and its child formulations

        Returns:
            bool: True if enough steps are initialized, False otherwise
        """
        return (self.GetBaseModelPart().ProcessInfo[Kratos.STEP] + 1 >=
                self.GetMinimumBufferSize())

    def IsConverged(self):
        """Recursively checks whether all formulations are converged.

        Returns:
            bool: True if all of them have converged, False if not
        """
        for formulation in self.__list_of_formulations:
            if (not formulation.IsConverged()):
                return False

        if (self.GetStrategy() is not None):
            is_converged = self.GetStrategy().IsConverged()
            if (is_converged):
                Kratos.Logger.PrintInfo(self.__class__.__name__, " *** CONVERGENCE ACHIEVED ***")
            return is_converged

        return True

    def ElementHasNodalProperties(self):
        """Recursively checks whether any one of the formulations require nodal properties.

        Returns:
            bool: True if one of them require nodal properties
        """

        for formulation in self.__list_of_formulations:
            if (formulation.ElementHasNodalProperties()):
                return True

        return False

    def GetElementNames(self):
        """Returns list of element names used in this formulation

        This method returns element names used in this formulation and
        child formulations

        Returns:
            list: a string list containing all the element names
        """
        element_names_list = []
        for formulation in self.__list_of_formulations:
            element_names_list.extend(formulation.GetElementNames())

        return element_names_list

    def GetConditionNames(self):
        """Returns list of condition names used in this formulation

        This method returns condition names used in this formulation and
        child formulations

        Returns:
            list: a string list containing all the condition names
        """

        condition_names_list = []
        for formulation in self.__list_of_formulations:
            condition_names_list.extend(formulation.GetConditionNames())

        return condition_names_list

    def GetMoveMeshFlag(self):
        """Returns move mesh flag

        Returns:
            bool: True if mesh move, False if not
        """
        return self.__move_mesh

    def SetCommunicator(self, communicator):
        """Sets the communicator for MPI use
        """
        self.__ExecuteRansFormulationMethods("SetCommunicator", [communicator])
        self.__communicator = communicator

    def GetCommunicator(self):
        """Get the communicator for MPI use

        Returns:
            Kratos.Communicator: Communicator used in the model part
        """
        if hasattr(self, "_RansFormulation__communicator"):
            return self.__communicator
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetCommunicator\" first before retrieving.")

    def IsPeriodic(self):
        """Checks whether current formulations are solved with periodic conditions

        Returns:
            bool: True if Periodic, False if not
        """
        if hasattr(self, "_RansFormulation__is_periodic"):
            return self.__is_periodic
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetIsPeriodic\" first before checking.")

    def SetIsPeriodic(self, value):
        """Sets periodicity recursively for all formulations

        Args:
            value (bool): True if formulations needs to be Periodic, False otherwise
        """
        self.__ExecuteRansFormulationMethods("SetIsPeriodic", [value])
        self.__is_periodic = value

    def SetTimeSchemeSettings(self, settings):
        """Sets time scheme settings recursively for all formulations

        Args:
            settings (Kratos.Parameters): Time scheme settings
        """
        self.__ExecuteRansFormulationMethods("SetTimeSchemeSettings", [settings])
        self.__time_scheme_settings = settings

    def GetTimeSchemeSettings(self):
        """Returns time scheme settings

        Returns:
            Kratos.Parameters: Time scheme settings used for formulations
        """
        if (hasattr(self, "_RansFormulation__time_scheme_settings")):
            return self.__time_scheme_settings
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetTimeSchemeSettings\" first before calling \"GetTimeSchemeSettings\".")

    def SetWallFunctionSettings(self):
        """Sets wall function settings
        """
        self.__ExecuteRansFormulationMethods("SetWallFunctionSettings")

    def GetBaseModelPart(self):
        """Returns base model part used in the formulation

        Returns:
            Kratos.ModelPart: Base model part used in the formulation
        """
        return self.__base_computing_model_part

    def GetComputingModelPart(self):
        """Returns computing model part used in the formulation

        Returns:
            Kratos.ModelPart: Computing model part used in the formulation
        """
        if (self.__base_computing_model_part.HasSubModelPart("fluid_computational_model_part")):
            return self.__base_computing_model_part.GetSubModelPart("fluid_computational_model_part")
        else:
            raise Exception("Computing model part \"fluid_computational_model_part\" not found as a submodel part in " + self.__base_computing_model_part.Name() + ".")

    def SetMaxCouplingIterations(self, max_iterations):
        """Sets max coupling iterations

        This is not done recursively because, there are some formulations which doesn't have coupling iterations.
        Each formulation needs to set this seperately if base class SolveCouplingStep is used.

        Args:
            max_iterations (int): Maximum number of coupling iterations to be done in the child formulations
        """
        self.__max_coupling_iterations = max_iterations

    def GetMaxCouplingIterations(self):
        """Returns maxmum number of coupling iterations used in this formulation

        Returns:
            int: Maximum number of coupling iterations
        """
        if (hasattr(self, "_RansFormulation__max_coupling_iterations")):
            return self.__max_coupling_iterations
        else:
            raise Exception(self.__class__.__name__ + " needs to use \"SetMaxCouplingIterations\" first before calling \"GetMaxCouplingIterations\".")

    def SetConstants(self, settings):
        """Recursively sets constants in the formulations

        Args:
            settings (Kratos.Parameters): Constants settings
        """
        self.__ExecuteRansFormulationMethods("SetConstants", [settings])

    def GetRansFormulationsList(self):
        """Returns list of formulations in this formulation

        Returns:
            List(RansFormulation): List of formulations in this formulation
        """
        return self.__list_of_formulations

    def GetSolvingVariables(self):
        """Returns list of variables being solved in this formulation

        Returns:
            List(RansFormulation): List of variables
        """
        variables = []
        for formulation in self.__list_of_formulations:
            variables.extend(formulation.GetSolvingVariables())

        return variables

    def GetProcessList(self):
        """Returns list of processes used in this formulation

        Returns:
            List(Kratos.RANSApplication.RansFormulationProcess): List of rans formulation processes in this formulation
        """
        return self.__list_of_processes

    def GetModelPart(self):
        """Returns the model part used for solving current formulation (if a strategy is used only.)

        Returns:
            Kratos.ModelPart: Model part used for solving current formulation
        """

        if (self.GetStrategy() is not None):
            return self.GetStrategy().GetModelPart()
        return None

    def DeleteModelPartsRecursively(self):
        for formulation in self.__list_of_formulations:
            formulation.DeleteModelPartsRecursively()

        if (self.GetModelPart() is not None):
            self.GetModelPart().GetModel().DeleteModelPart(self.GetModelPart().FullName())

    def GetStrategy(self):
        """Returns strategy used in this formulation, if used any.

        Returns:
            Kratos.ImplicitSolvingStrategy: Strategy used in this formulation, None if not used.
        """
        return None

    def GetInfo(self):
        """Recursively identify formulations being used.

        Returns:
            str: Information of all the formulations
        """
        info = "\n" + self.__class__.__name__
        if (self.GetModelPart() is not None):
            info += "\n   Model part    : " + str(self.GetModelPart().FullName())

        if (self.GetMaxCouplingIterations() != 0):
            info += "\n   Max iterations: " + str(self.GetMaxCouplingIterations())

        if (len(self.GetProcessList()) > 0):
            info += "\n   Process list:"
            for process in self.GetProcessList():
                info += "\n      " + str(process).strip()

        for formulation in self.GetRansFormulationsList():
            info += str(formulation.GetInfo()).replace("\n", "\n   ")
        return info

    def __ExecuteRansFormulationMethods(self, method_name, args = []):
        for formulation in self.__list_of_formulations:
            getattr(formulation, method_name)(*args)

    def __ExecuteProcessMethods(self, method_name):
        for process in self.__list_of_processes:
            getattr(process, method_name)()




