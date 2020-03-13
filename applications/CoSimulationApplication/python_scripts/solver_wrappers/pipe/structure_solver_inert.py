import numpy as np
import os.path as path
from scipy.linalg import solve_banded

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperPipeStructure(parameters)


class SolverWrapperPipeStructure(CoSimulationComponent):
    Al = 2  # Number of terms below diagonal in matrix
    Au = 2  # Number of terms above diagonal in matrix

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
        d = self.settings["d"].GetDouble()  # Diameter
        self.r0 = d / 2.0  # Reference radius of cross section
        self.rhof = self.settings["rhof"].GetDouble()  # Fluid density

        self.preference = self.settings["preference"].GetDouble() if self.settings.Has("preference") else 0.0  # Reference pressure

        e = self.settings["e"].GetDouble()  # Young's modulus of structure
        nu = self.settings["nu"].GetDouble()  # Poisson's ratio
        self.h = self.settings["h"].GetDouble()  # Thickness of structure
        self.rhos = self.settings["rhos"].GetDouble()  # Structure density
        self.cmk2 = (e * self.h) / (self.rhof * d)  # Wave speed squared
        self.b1 = (self.h * e) / (1 - nu ** 2) * (self.h ** 2) / 12
        self.b2 = self.b1 * (2 * nu) / self.r0 ** 2
        self.b3 = (self.h * e) / (1 - nu ** 2) * 1 / self.r0 ** 2

        self.m = self.settings["m"].GetInt()  # Number of segments
        self.dz = l / self.m  # Segment length
        self.z = np.arange(self.dz / 2.0, l, self.dz)  # Data is stored in cell centers

        self.k = 0  # Iteration
        self.n = 0  # Time step (no restart implemented)
        self.dt = self.settings["delta_t"].GetDouble()  # Time step size

        self.gamma = self.settings["gamma"].GetDouble()  # Newmark parameter: gamma >= 1/2
        self.beta = self.settings["beta"].GetDouble()  # Newmark parameter: beta >= 1/4 * (1/2 + gamma) ^ 2
        if not self.gamma >= 0.5 or not self.beta >= 0.25 * (0.5 + self.gamma) ** 2:
            raise Exception("Inadequate Newmark parameteres")

        self.newtonmax = self.settings["newtonmax"].GetInt()  # Maximal number of Newton iterations
        self.newtontol = self.settings["newtontol"].GetDouble()  # Tolerance of Newton iterations

        # Initialization
        self.areference = np.pi * self.r0 ** 2  # Reference area of cross section
        self.p = np.ones(self.m) * self.preference  # Pressure
        self.a = np.ones(self.m) * self.areference  # Area of cross section
        self.r = np.ones(self.m + 4) * self.r0  # Radius of cross section
        self.rn = np.ones(self.m + 4) * self.r0  # Previous radius of cross section

        self.rdot = np.zeros(self.m)
        self.rddot = np.zeros(self.m)
        self.rndot = np.zeros(self.m)
        self.rnddot = np.zeros(self.m)

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
            node.SetSolutionStepValue(self.variable_pres, step, self.p[0])
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
        self.rn = np.array(self.r)
        self.rndot = np.array(self.rdot)
        self.rnddot = np.array(self.rddot)

    def SolveSolutionStep(self, interface_input):
        # Input does not contain boundary conditions
        p = interface_input.GetNumpyArray()
        self.interface_input.SetNumpyArray(p)
        self.p = np.array(p)

        # Newton iterations
        converged = False
        f = self.GetResidual()
        residual0 = np.linalg.norm(f)
        if residual0:
            for s in range(self.newtonmax):
                j = self.GetJacobian()
                b = -f
                x = solve_banded((self.Al, self.Au), j, b)
                self.r += x
                f = self.GetResidual()
                residual = np.linalg.norm(f)
                if residual / residual0 < self.newtontol:
                    converged = True
                    break
            if not converged:
                Exception("Newton failed to converge")

        self.k += 1
        if self.debug:
            file_name = self.working_directory + f"/Area_TS{self.n}_IT{self.k}"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

        # Output does not contain boundary conditions
        self.a = self.r[2:self.m + 2] ** 2 * np.pi
        self.interface_output.SetNumpyArray(self.a)
        return self.interface_output.deepcopy()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.rddot = ((self.r[2:self.m + 2] - self.rn[2:self.m + 2]) / (self.beta * self.dt ** 2)
                      - self.rndot / (self.beta * self.dt) - self.rnddot * (1 / (2 * self.beta) - 1))
        self.rdot = self.rndot + self.dt * (1 - self.gamma) * self.rnddot + self.dt * self.gamma * self.rddot

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

    def GetResidual(self):
        f = np.zeros(self.m + 4)
        f[0] = self.r[0] - self.r0
        f[1] = self.r[1] - self.r0
        f[2:self.m + 2] = ((self.rhos * self.h) / (self.beta * self.dt ** 2) * self.r[2: self.m + 2]
                           + self.b1 / self.dz ** 4 * (self.r[4:self.m + 4] - 4.0 * self.r[3:self.m + 3]
                                                       + 6.0 * self.r[2:self.m + 2] - 4.0 * self.r[1:self.m + 1]
                                                       + self.r[0:self.m])
                           - self.b2 / self.dz ** 2 * (self.r[3:self.m + 3] - 2.0 * self.r[2:self.m + 2]
                                                       + self.r[1:self.m + 1])
                           + self.b3 * (self.r[2:self.m + 2] - self.r0)
                           - (self.p - self.preference)
                           - self.rhos * self.h * (self.rn[2:self.m + 2] / (self.beta * self.dt ** 2)
                                                   + self.rndot / (self.beta * self.dt)
                                                   + self.rnddot * (1.0 / (2.0 * self.beta) - 1.0)))
        f[self.m + 2] = self.r[self.m + 2] - self.r0
        f[self.m + 3] = self.r[self.m + 3] - self.r0
        return f

    def GetJacobian(self):
        j = np.zeros((self.Al + self.Au + 1, self.m + 4))
        j[self.Au + 0 - 0, 0] = 1.0  # [0, 0]
        j[self.Au + 1 - 1, 1] = 1.0  # [1, 1]
        j[self.Au + 2, 0:self.m] = self.b1 / self.dz ** 4  # [i, (i - 2)]
        j[self.Au + 1, 1:self.m + 1] = - 4.0 * self.b1 / self.dz ** 4 - self.b2 / self.dz ** 2  # [i, (i - 1)]
        j[self.Au + 0, 2:self.m + 2] = ((self.rhos * self.h) / (self.beta * self.dt ** 2)
                                        + 6.0 * self.b1 / self.dz ** 4 + 2.0 * self.b2 / self.dz ** 2
                                        + self.b3)  # [i, i]
        j[self.Au - 1, 3:self.m + 3] = - 4.0 * self.b1 / self.dz ** 4 - self.b2 / self.dz ** 2  # [i, (i + 1)]
        j[self.Au - 2, 4:self.m + 4] = self.b1 / self.dz ** 4  # [i, (i + 2)]
        j[self.Au + (self.m + 2) - (self.m + 2), self.m + 2] = 1.0  # [m + 2, m + 2]
        j[self.Au + (self.m + 3) - (self.m + 3), self.m + 3] = 1.0  # [m + 3, m + 3]
        return j
