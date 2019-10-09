import numpy as np
import math as m
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
        working_directory = self.settings["working_directory"].GetString()
        input_file = self.settings["input_file"].GetString()
        settings_file_name = path.join(working_directory, input_file)
        with open(settings_file_name, 'r') as settings_file:
            self.settings.AddMissingParameters(cs_data_structure.Parameters(settings_file.read()))

        # Settings
        l = self.settings["l"].GetDouble()  # Length
        self.d = self.settings["d"].GetDouble() # Diameter
        self.rhof = self.settings["rhof"].GetDouble()  # Density

        e = self.settings["e"].GetDouble()  # Young"s modulus of structure
        h = self.settings["h"].GetDouble()  # Thickness of structure
        self.cmk2 = (e * h) / (self.rhof * self.d)  # Wave speed squared

        self.m = self.settings["m"].GetInt()  # Number of segments
        self.dz = l / self.m  # Segment length
        self.z = np.arange(self.dz / 2.0, l, self.dz)  # Data is stored in cell centers

        self.n = 0  # Time step

        # Initialization
        self.p = np.ones(self.m) * 2.0 * self.cmk2  # Pressure
        self.a = np.ones(self.m) * m.pi * self.d ** 2 / 4.0  # Area of cross section
        self.p0 = 0.0  # Reference pressure
        self.a0 = m.pi * self.d ** 2 / 4.0  # Reference area of cross section
        self.c02 = self.cmk2 - self.p0 / 2.0  # Wave speed squared with reference pressure

        # ModelParts
        self.variable_pres = KM.KratosGlobals.GetVariable("PRESSURE")
        self.variable_area = KM.KratosGlobals.GetVariable("AREA")
        self.model = cs_data_structure.Model()
        self.model_part = self.model.CreateModelPart("wall")
        self.model_part.AddNodalSolutionStepVariable(self.variable_pres)
        self.model_part.AddNodalSolutionStepVariable(self.variable_area)
        for i in range(len(self.z)):
            self.model_part.CreateNewNode(i, 0.0, 0.0, self.z[i])
        step = 0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(self.variable_pres, step, self.p[0])
            node.SetSolutionStepValue(self.variable_area, step, self.a[0])

        # Interfaces
        self.interface_input = CoSimulationInterface(self.model, self.settings["interface_input"])
        self.interface_output = CoSimulationInterface(self.model, self.settings["interface_output"])

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.n += 1

    def SolveSolutionStep(self, interface_input):
        self.interface_input = interface_input
        self.p = self.interface_input.GetNumpyArray()

        # Independent rings model
        for i in range(len(self.p)):
            if self.p[i] > 2.0 * self.c02 + self.p0:
                raise ValueError("Unphysical pressure")
        for i in range(len(self.a)):
            self.a[i] = self.a0 * (2.0 / (2.0 + (self.p0 - self.p[i]) / self.c02)) ** 2

        self.interface_output.SetNumpyArray(self.a)
        return self.interface_output

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

    def GetInterfaceInput(self):
        return self.interface_input

    def SetInterfaceInput(self):
        Exception("This solver wrapper provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output

    def SetInterfaceOutput(self):
        Exception("This solver wrapper provides no mapping.")
