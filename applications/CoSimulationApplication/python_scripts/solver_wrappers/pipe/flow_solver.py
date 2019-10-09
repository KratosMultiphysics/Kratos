import numpy as np
import math as m
import os.path as path
from scipy.linalg import solve_banded

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperPipeFlow(parameters)


class SolverWrapperPipeFlow(CoSimulationComponent):
    Al = 4  # Number of terms below diagonal in matrix
    Au = 4  # Number of terms above diagonal in matrix

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
        self.d = self.settings["d"].GetDouble()  # Diameter
        self.rhof = self.settings["rhof"].GetDouble()  # Density

        self.ureference = self.settings["ureference"].GetDouble()  # Reference of inlet boundary condition
        self.uamplitude = self.settings["uamplitude"].GetDouble()  # Amplitude of inlet boundary condition
        self.uperiod = self.settings["uperiod"].GetDouble()  # Period of inlet boundary condition
        self.utype = self.settings["utype"].GetInt()  # Type of inlet boundary condition

        e = self.settings["e"].GetDouble()  # Young"s modulus of structure
        h = self.settings["h"].GetDouble()  # Thickness of structure
        self.cmk2 = (e * h) / (self.rhof * self.d)  # Wave speed squared of outlet boundary condition

        self.m = self.settings["m"].GetInt()  # Number of segments
        self.dz = l / self.m  # Segment length
        self.z = np.arange(self.dz / 2.0, l, self.dz)  # Data is stored in cell centers

        self.n = 0  # Time step
        self.dt = self.settings["dt"].GetDouble()  # Time step size
        self.alpha = 0.0  # Numerical damping parameter due to central discretization of pressure in momentum equation

        self.newtonmax = self.settings["newtonmax"].GetInt()  # Maximal number of Newton iterations
        self.newtontol = self.settings["newtontol"].GetDouble()  # Tolerance of Newton iterations

        # Initialization
        self.u = np.ones(self.m + 2) * self.ureference  # Velocity
        self.un = np.ones(self.m + 2) * self.ureference  # Previous velocity
        self.p = np.zeros(self.m + 2)  # Pressure
        self.pn = np.zeros(self.m + 2)  # Previous pressure (only value at outlet is used)
        self.a = np.ones(self.m + 2) * m.pi * self.d ** 2 / 4.0  # Area of cross section
        self.an = np.ones(self.m + 2) * m.pi * self.d ** 2 / 4.0  # Previous area of cross section

        # ModelParts
        self.variable_area = KM.KratosGlobals.GetVariable("AREA")
        self.variable_pres = KM.KratosGlobals.GetVariable("PRESSURE")
        self.model = cs_data_structure.Model()
        self.model_part = self.model.CreateModelPart("wall")
        self.model_part.AddNodalSolutionStepVariable(self.variable_area)
        self.model_part.AddNodalSolutionStepVariable(self.variable_pres)
        for i in range(len(self.z)):
            self.model_part.CreateNewNode(i, 0.0, 0.0, self.z[i])
        step = 0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(self.variable_area, step, self.a[0])
            node.SetSolutionStepValue(self.variable_pres, step, self.p[0])

        # Interfaces
        self.interface_input = CoSimulationInterface(self.model, self.settings["interface_input"])
        self.interface_output = CoSimulationInterface(self.model, self.settings["interface_output"])

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.n += 1
        self.un = np.array(self.u)
        self.pn = np.array(self.p)
        self.an = np.array(self.a)

    def SolveSolutionStep(self, interface_input):
        # Input does not contain boundary conditions
        self.interface_input = interface_input
        a = self.interface_input.GetNumpyArray()
        self.a[1:self.m + 1] = a
        self.a[0] = self.a[1]
        self.a[self.m + 1] = self.a[self.m]

        # Newton iterations
        converged = False
        f = self.GetResidual()
        residual0 = np.linalg.norm(f)
        if residual0:
            for s in range(self.newtonmax):
                j = self.GetJacobian()
                b = -f
                x = solve_banded((self.Al, self.Au), j, b)
                self.u += x[0::2]
                self.p += x[1::2]
                self.u[0] = self.GetBoundary()
                f = self.GetResidual()
                residual = np.linalg.norm(f)
                if residual / residual0 < self.newtontol:
                    converged = True
                    break
            if not converged:
                Exception("Newton failed to converge")

        # Output does not contain boundary conditions
        p = self.p[1:self.m + 1]
        self.interface_output.SetNumpyArray(p)
        return self.interface_output

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

    def GetInterfaceInput(self):
        return self.interface_input

    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output

    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    def GetBoundary(self):
        if self.utype == 1:
            u = self.ureference + self.uamplitude * m.sin(2.0 * m.pi * (self.n * self.dt) / self.uperiod)
        elif self.utype == 2:
            u = self.ureference + self.uamplitude
        elif self.utype == 3:
            u = self.ureference + self.uamplitude * (m.sin(m.pi * (self.n * self.dt) / self.uperiod)) ** 2
        else:
            u = self.ureference + self.uamplitude * (self.n * self.dt) / self.uperiod
        return u

    def GetResidual(self):
        usign = self.u[1:self.m + 1] > 0
        ur = self.u[1:self.m + 1] * usign + self.u[2:self.m + 2] * (1.0 - usign)
        ul = self.u[0:self.m] * usign + self.u[1:self.m + 1] * (1.0 - usign)
        self.alpha = m.pi * self.d ** 2 / 4.0 / (self.ureference + self.dz / self.dt)

        f = np.zeros(2 * self.m + 4)
        f[0] = self.u[0] - self.GetBoundary()
        f[1] = self.p[0] - (2.0 * self.p[1] - self.p[2])
        f[2:2 * self.m + 2:2] = (self.dz / self.dt * (self.a[1:self.m + 1] - self.an[1:self.m + 1])
                                 + (self.u[1:self.m + 1] + self.u[2:self.m + 2])
                                 * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                 - (self.u[1:self.m + 1] + self.u[0:self.m])
                                 * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0
                                 - self.alpha * (self.p[2:self.m + 2] - 2.0 * self.p[1:self.m + 1] + self.p[0:self.m]))
        f[3:2 * self.m + 3:2] = (self.dz / self.dt * (self.u[1:self.m + 1] * self.a[1:self.m + 1]
                                 - self.un[1:self.m + 1] * self.an[1:self.m + 1])
                                 + ur * (self.u[1:self.m + 1] + self.u[2:self.m + 2])
                                 * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                 - ul * (self.u[1:self.m + 1] + self.u[0:self.m])
                                 * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0
                                 + ((self.p[2:self.m + 2] - self.p[1:self.m + 1])
                                 * (self.a[1:self.m + 1] + self.a[2:self.m + 2])
                                 + (self.p[1:self.m + 1] - self.p[0:self.m])
                                 * (self.a[1:self.m + 1] + self.a[0:self.m])) / 4.0)
        f[2 * self.m + 2] = self.u[self.m + 1] - (2.0 * self.u[self.m] - self.u[self.m - 1])
        f[2 * self.m + 3] = self.p[self.m + 1] - (2.0 * (self.cmk2 - (m.sqrt(self.cmk2 - self.pn[self.m + 1] / 2.0)
                             - (self.u[self.m + 1] - self.un[self.m + 1])/ 4.0) ** 2))
        return f

    def GetJacobian(self):
        usign = self.u[1:self.m + 1] > 0
        j = np.zeros((self.Al + self.Au + 1, 2 * self.m + 4))
        j[self.Au + 0 - 0, 0] = 1.0  # [0,0]
        j[self.Au + 1 - 1, 1] = 1.0  # [1,1]
        j[self.Au + 1 - 3, 3] = -2.0  # [1,3]
        j[self.Au + 1 - 5, 5] = 1.0  # [1,5]

        j[self.Au + 2, 0:2 * self.m + 0:2] = -(self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0  # [2*i, 2*(i-1)]
        j[self.Au + 3, 0:2 * self.m + 0:2] = (-((self.u[1:self.m + 1] + 2.0 * self.u[0:self.m]) * usign
                                                + self.u[1:self.m + 1] * (1.0 - usign))
                                                * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0)  # [2*i+1, 2*(i-1)]
        j[self.Au + 1, 1:2 * self.m + 1:2] = -self.alpha  # [2*i, 2*(i-1)+1]
        j[self.Au + 2, 1:2 * self.m + 1:2] = -(self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0  # [2*i+1, 2*(i-1)+1]

        j[self.Au + 0, 2:2 * self.m + 2:2] = ((self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                                - (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0)  # [2*i, 2*i]
        j[self.Au + 1, 2:2 * self.m + 2:2] = (self.dz / self.dt * self.a[1:self.m + 1]
                                                  + ((2.0 * self.u[1:self.m + 1] + self.u[2:self.m + 2]) * usign
                                                     + self.u[2:self.m + 2] * (1.0 - usign))
                                                  * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                                  - (self.u[0:self.m] * usign
                                                     + (2.0 * self.u[1:self.m + 1] + self.u[0:self.m]) * (1.0 - usign))
                                                  * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0)  # [2*i+1, 2*i]
        j[self.Au - 1, 3:2 * self.m + 3:2] = 2.0 * self.alpha  # [2*i, 2*i+1]
        j[self.Au + 0, 3:2 * self.m + 3:2] = (-(self.a[1:self.m + 1] + self.a[2:self.m + 2])
                                                  + (self.a[1:self.m + 1] + self.a[0:self.m])) / 4.0  # [2*i+1, 2*i+1]

        j[self.Au - 2, 4:2 * self.m + 4:2] = (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0  # [2*i, 2*(i+1)]
        j[self.Au - 1, 4:2 * self.m + 4:2] = ((self.u[1:self.m + 1] * usign + (self.u[1:self.m + 1]
                                                                                   + 2.0 * self.u[2:self.m + 2])
                                                   * (1.0 - usign))
                                                  * (self.a[1:self.m + 1]
                                                     + self.a[2:self.m + 2]) / 4.0)  # [2*i+1, 2*(i+1)]
        j[self.Au - 3, 5:2 * self.m + 5:2] = -self.alpha  # [2*i, 2*(i+1)+1]
        j[self.Au - 2, 5:2 * self.m + 5:2] = (self.a[1:self.m + 1]
                                                  + self.a[2:self.m + 2]) / 4.0  # [2*i+1, 2*(i+1)+1]

        j[self.Au + (2 * self.m + 2) - (2 * self.m + 2), 2 * self.m + 2] = 1.0  # [2*m+2, 2*m+2]
        j[self.Au + (2 * self.m + 2) - (2 * self.m), 2 * self.m] = -2.0  # [2*m+2, 2*m]
        j[self.Au + (2 * self.m + 2) - (2 * self.m - 2), 2 * self.m - 2] = 1.0  # [2*m+2, 2*m-2]
        j[self.Au + (2 * self.m + 3) - (2 * self.m + 2), 2 * self.m + 2] = (-(m.sqrt(self.cmk2
                                                                                         - self.pn[self.m + 1] / 2.0)
                                                                                  - (self.u[self.m + 1]
                                                                                     - self.un[self.n + 1])
                                                                                  / 4.0))  # [2*m+3, 2*m+2]
        j[self.Au + (2 * self.m + 3) - (2 * self.m + 3), 2 * self.m + 3] = 1.0  # [2*m+3, 2*m+3]

        return j
