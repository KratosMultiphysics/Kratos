import numpy as np
import math as m

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverInterfacePipeStructure(parameters)


class SolverInterfacePipeStructure(CoSimulationComponent):
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
        self.dt = 0  # Time step size

        # Initialization
        self.p = np.ones(self.m) * 2.0 * self.cmk2  # Pressure
        self.a = np.ones(self.m) * m.pi * self.d ** 2 / 4.0  # Area of cross section
        self.p0 = 0.0  # Reference pressure
        self.a0 = m.pi * self.d ** 2 / 4.0  # Reference area of cross section
        self.c02 = self.cmk2 - self.p0 / 2.0  # Wave speed squared with reference pressure

        # ModelPart of interface
        self.model = cs_data_structure.Model()
        self.model_part = self.model = self.model.CreateModelPart("pipe_flow")
        self.variable_pres = cs_data_structure.VariableDouble("PRESSURE")
        self.model_part.AddNodalSolutionStepVariable(self.variable_pres)
        self.variable_disp = cs_data_structure.VariableDouble("DISPLACEMENT_Y")
        self.model_part.AddNodalSolutionStepVariable(self.variable_disp)
        for i in range(len(self.z)):
            self.model_part.CreateNewNode(i, 0.0, 0.0, self.z[i])

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.n += 1

    def calculate(self, p):
        # Independent rings model
        self.p = p
        for i in range(len(self.p)):
            if self.p[i] > 2.0 * self.c02 + self.p0:
                raise ValueError("Unphysical pressure")
        for i in range(len(self.a)):
            self.a[i] = self.a0 * (2.0 / (2.0 + (self.p0 - self.p[i]) / self.c02)) ** 2
        # Return copy of output
        return np.array(self.a)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()
