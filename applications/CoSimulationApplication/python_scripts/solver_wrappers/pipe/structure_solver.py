import numpy as np
import math as m

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperPipeStructure(parameters)


class SolverWrapperPipeStructure(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

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
        self.model = cs_data_structure.Model()
        self.model_part_in = self.model.CreateModelPart("pipe_structure_in")
        self.variable_pres = cs_data_structure.VariableDouble("PRESSURE")
        self.model_part_in.AddNodalSolutionStepVariable(self.variable_pres)
        self.model_part_out = self.model.CreateModelPart("pipe_structure_out")
        self.variable_area = cs_data_structure.VariableDouble("AREA")
        self.model_part_out.AddNodalSolutionStepVariable(self.variable_area)
        for i in range(len(self.z)):
            self.model_part_in.CreateNewNode(i, 0.0, 0.0, self.z[i])
            self.model_part_out.CreateNewNode(i, 0.0, 0.0, self.z[i])
        step = 0
        for node in self.model_part_in.Nodes:
            node.SetSolutionStepValue(self.variable_pres, step, self.p[0])
        for node in self.model_part_out.Nodes:
            node.SetSolutionStepValue(self.variable_area, step, self.a[0])

        # Interfaces
        self.interface_in = CoSimulationInterface([self.model_part_in])
        self.interface_out = CoSimulationInterface([self.model_part_out])

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.n += 1

    def SolveSolutionStep(self, interface_in):
        self.interface_in = interface_in
        self.p = self.interface_in.GetNumpyArray()

        # Independent rings model
        for i in range(len(self.p)):
            if self.p[i] > 2.0 * self.c02 + self.p0:
                raise ValueError("Unphysical pressure")
        for i in range(len(self.a)):
            self.a[i] = self.a0 * (2.0 / (2.0 + (self.p0 - self.p[i]) / self.c02)) ** 2

        self.interface_out.SetNumpyArray(self.a)
        return self.interface_out

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

    def GetInterfaceIn(self):
        return self.interface_in

    def SetInterfaceIn(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOut(self):
        return self.interface_out

    def SetInterfaceOut(self):
        Exception("This solver interface provides no mapping.")
