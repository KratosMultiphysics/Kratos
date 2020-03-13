import numpy as np
import os.path as path

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperPipeStructure(parameters)


class SolverWrapperPipeStructure(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        # Reading
        self.parameters = parameters
        self.settings = parameters["settings"]
        self.working_directory = self.settings["working_directory"].GetString()
        input_file = self.settings["input_file"].GetString()
        settings_file_name = path.join(self.working_directory, input_file)
        with open(settings_file_name, 'r') as settings_file:
            self.settings.AddMissingParameters(cs_data_structure.Parameters(settings_file.read()))

        # Settings
        l = self.settings["l"].GetDouble()  # Length
        self.d = self.settings["d"].GetDouble()  # Diameter
        self.rhof = self.settings["rhof"].GetDouble()  # Fluid density

        self.preference = self.settings["preference"].GetDouble() if self.settings.Has("preference") else 0.0  # Reference pressure

        e = self.settings["e"].GetDouble()  # Young's modulus of structure
        h = self.settings["h"].GetDouble()  # Thickness of structure
        self.cmk2 = (e * h) / (self.rhof * self.d)  # Wave speed squared

        self.m = self.settings["m"].GetInt()  # Number of segments
        self.dz = l / self.m  # Segment length
        self.z = np.arange(self.dz / 2.0, l, self.dz)  # Data is stored in cell centers

        self.k = 0  # Iteration
        self.n = 0  # Time step (no restart implemented)

        # Initialization
        self.areference = np.pi * self.d ** 2 / 4  # Reference area of cross section
        self.p = np.ones(self.m) * self.preference  # Kinematic pressure
        self.a = np.ones(self.m) * self.areference  # Area of cross section
        self.c02 = self.cmk2 - self.preference / 2.0  # Wave speed squared with reference pressure

        # ModelParts
        self.variable_pres = vars(KM)["PRESSURE"]
        self.variable_area = vars(KM)["AREA"]
        self.model = cs_data_structure.Model()
        self.model_part = self.model.CreateModelPart("wall")
        self.model_part.AddNodalSolutionStepVariable(self.variable_pres)
        self.model_part.AddNodalSolutionStepVariable(self.variable_area)
        for i in range(len(self.z)):
            self.model_part.CreateNewNode(i, 0.0, 0.0, self.z[i])
        step = 0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(self.variable_pres, step, self.p[0] * self.rhof)
            node.SetSolutionStepValue(self.variable_area, step, self.a[0])

        # Interfaces
        self.interface_input = CoSimulationInterface(self.model, self.settings["interface_input"])
        self.interface_output = CoSimulationInterface(self.model, self.settings["interface_output"])

        # Debug
        self.debug = True  # Set on true to save solution of each time step
        self.OutputSolutionStep()

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.k = 0
        self.n += 1

    def SolveSolutionStep(self, interface_input):
        self.interface_input = interface_input.deepcopy()
        self.p = self.interface_input.GetNumpyArray() / self.rhof

        # Independent rings model
        for i in range(len(self.p)):
            if self.p[i] > 2.0 * self.c02 + self.preference:
                raise ValueError("Unphysical pressure")
        self.a = self.areference * (2.0 / (2.0 + (self.preference - self.p) / self.c02)) ** 2

        self.k += 1
        if self.debug:
            file_name = self.working_directory + f"/Area_TS{self.n}_IT{self.k}"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

        self.interface_output.SetNumpyArray(self.a)
        return self.interface_output.deepcopy()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

    def OutputSolutionStep(self):
        if self.debug:
            file_name = self.working_directory + f"/Area_TS{self.n}"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()

    def SetInterfaceInput(self):
        raise Exception("This solver wrapper provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()

    def SetInterfaceOutput(self):
        raise Exception("This solver wrapper provides no mapping.")
