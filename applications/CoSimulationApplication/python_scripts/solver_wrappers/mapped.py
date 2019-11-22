from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperMapped(parameters)


class SolverWrapperMapped(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        # Read parameters
        self.parameters = parameters
        self.settings = parameters["settings"]

        # Create solver
        self.solver_wrapper = cs_tools.CreateInstance(self.settings["solver_wrapper"])

    def Initialize(self):
        super().Initialize()

        self.solver_wrapper.Initialize()

        self.mapper_interface_input = None
        self.mapper_interface_output = None

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.solver_wrapper.InitializeSolutionStep()

    def SolveSolutionStep(self, interface_input_from):
        if self.mapper_interface_input is None:
            self.SetInterfaceInput(interface_input_from)
        if self.mapper_interface_output is None:
            self.SetInterfaceOutput(interface_input_from)

        self.interface_input_from = interface_input_from
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        self.interface_output_from = self.solver_wrapper.SolveSolutionStep(self.interface_input_to)
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.solver_wrapper.FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

        self.solver_wrapper.Finalize()
        self.mapper_interface_input.Finalize()
        self.mapper_interface_output.Finalize()

    def GetInterfaceInput(self):
        # Does not contain most recent data
        return self.interface_input_from

    def SetInterfaceInput(self, interface_input_from):
        # Create input mapper
        self.interface_input_from = interface_input_from
        self.interface_input_to = self.solver_wrapper.GetInterfaceInput()
        self.mapper_interface_input = cs_tools.CreateInstance(self.settings["mapper_interface_input"])
        self.mapper_interface_input.Initialize(self.interface_input_from, self.interface_input_to)

    def GetInterfaceOutput(self):
        return self.interface_output_to

    def SetInterfaceOutput(self, interface_output_to):
        # Create output mapper
        self.interface_output_to = interface_output_to
        self.interface_output_from = self.solver_wrapper.GetInterfaceOutput()
        self.mapper_interface_output = cs_tools.CreateInstance(self.settings["mapper_interface_output"])
        self.mapper_interface_output.Initialize(self.interface_output_from, self.interface_output_to)

