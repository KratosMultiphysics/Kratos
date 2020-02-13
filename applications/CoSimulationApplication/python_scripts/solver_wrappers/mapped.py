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

        solver_wrapper_parameters = self.settings["solver_wrapper"]
        solver_wrapper_settings = solver_wrapper_parameters["settings"]
        # Add time_step_start and delta_t to solver_wrapper_settings
        solver_wrapper_settings.AddValue("timestep_start", self.settings["timestep_start"])
        solver_wrapper_settings.AddValue("delta_t", self.settings["delta_t"])
        solver_wrapper_parameters.RemoveValue("settings")
        solver_wrapper_parameters.AddValue("settings", solver_wrapper_settings)

        # Create solver
        self.solver_wrapper = cs_tools.CreateInstance(self.settings["solver_wrapper"])

    def Initialize(self):
        super().Initialize()

        self.solver_wrapper.Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.solver_wrapper.InitializeSolutionStep()

    def SolveSolutionStep(self, interface_input_from):
        self.interface_input_from = interface_input_from.deepcopy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        self.interface_output_from = self.solver_wrapper.SolveSolutionStep(self.interface_input_to)
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to.deepcopy()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.solver_wrapper.FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

        self.solver_wrapper.Finalize()
        self.mapper_interface_input.Finalize()
        self.mapper_interface_output.Finalize()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        self.solver_wrapper.OutputSolutionStep()
        self.mapper_interface_input.OutputSolutionStep()
        self.mapper_interface_output.OutputSolutionStep()

    def GetInterfaceInput(self):
        # Does not contain most recent data
        return self.interface_input_from.deepcopy()

    def SetInterfaceInput(self, interface_input_from):
        # Create input mapper
        self.interface_input_from = interface_input_from.deepcopy()
        self.interface_input_to = self.solver_wrapper.GetInterfaceInput()

        self.mapper_interface_input = cs_tools.CreateInstance(self.settings["mapper_interface_input"])
        self.mapper_interface_input.Initialize(self.interface_input_from, self.interface_input_to)

    def GetInterfaceOutput(self):
        self.interface_output_from = self.solver_wrapper.GetInterfaceOutput()
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to.deepcopy()

    def SetInterfaceOutput(self, interface_output_to):
        # Create output mapper
        self.interface_output_to = interface_output_to.deepcopy()
        self.interface_output_from = self.solver_wrapper.GetInterfaceOutput()

        self.mapper_interface_output = cs_tools.CreateInstance(self.settings["mapper_interface_output"])
        self.mapper_interface_output.Initialize(self.interface_output_from, self.interface_output_to)

    def PrintInfo(self, label):
        cs_tools.Print(label, "The component ", self.__class__.__name__, "  maps the following solver wrapper:")
        self.solver_wrapper.PrintInfo(label + '\t')
        label += '\t'
        cs_tools.Print(label, "With input mapper:")
        self.mapper_interface_input.PrintInfo(label + '\t')
        cs_tools.Print(label, "And output mapper:")
        self.mapper_interface_output.PrintInfo(label + '\t')
