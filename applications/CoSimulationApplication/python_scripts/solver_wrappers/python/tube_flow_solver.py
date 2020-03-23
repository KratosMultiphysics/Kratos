import numpy as np
import os.path as path
from scipy.linalg import solve_banded

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperTubeFlow(parameters)


class SolverWrapperTubeFlow(CoSimulationComponent):
    Al = 4  # Number of terms below diagonal in matrix
    Au = 4  # Number of terms above diagonal in matrix

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
        self.rhof = self.settings["rhof"].GetDouble()  # Density

        self.ureference = self.settings["ureference"].GetDouble()  # Reference velocity
        self.preference = self.settings["preference"].GetDouble() if self.settings.Has("preference") else 0.0 # Reference pressure
        self.inlet_boundary = self.settings["inlet_boundary"]
        self.inlet_variable = self.inlet_boundary["variable"].GetString()  # Variable upon which boundary condition is specified
        if self.inlet_variable == "velocity":
            self.inlet_reference = self.inlet_boundary["reference"].GetDouble() if self.inlet_boundary.Has(
                "reference") else self.ureference  # Reference of velocity inlet boundary condition
        elif self.inlet_variable == "pressure":
            self.inlet_reference = self.inlet_boundary["reference"].GetDouble() if self.inlet_boundary.Has(
                "reference") else self.preference  # Reference of pressure inlet boundary condition
        else:
            raise ValueError(f"The inlet_variable \'{self.inlet_variable}\' is not implemented,"
                             f" chose between \'pressure\' and \'velocity\'")
        self.inlet_type = self.inlet_boundary["type"].GetInt()  # Type of inlet boundary condition
        self.inlet_amplitude = self.inlet_boundary["amplitude"].GetDouble()  # Amplitude of inlet boundary condition
        self.inlet_period = self.inlet_boundary["period"].GetDouble()  # Period of inlet boundary condition
        # Adjust to kinematic pressure
        if self.inlet_variable == "pressure":
            self.inlet_reference = self.inlet_reference / self.rhof
            self.inlet_amplitude = self.inlet_amplitude / self.rhof
        self.preference = self.preference / self.rhof

        self.outlet_boundary = self.settings["outlet_boundary"]
        self.outlet_type = self.outlet_boundary["type"].GetInt()  # Type of outlet boundary condition

        e = self.settings["e"].GetDouble()  # Young"s modulus of structure
        h = self.settings["h"].GetDouble()  # Thickness of structure
        self.cmk2 = (e * h) / (self.rhof * self.d)  # Wave speed squared of outlet boundary condition

        self.m = self.settings["m"].GetInt()  # Number of segments
        self.dz = l / self.m  # Segment length
        axial_offset = self.settings["axial_offset"].GetDouble() if self.settings.Has("axial_offset") else 0.0  # Start position along axis
        self.z = axial_offset + np.arange(self.dz / 2.0, l, self.dz)  # Data is stored in cell centers

        self.k = 0  # Iteration
        self.n = 0  # Time step (no restart implemented)
        self.dt = self.settings["delta_t"].GetDouble()  # Time step size
        self.alpha = 0.0  # Numerical damping parameter due to central discretization of pressure in momentum equation

        self.newtonmax = self.settings["newtonmax"].GetInt()  # Maximal number of Newton iterations
        self.newtontol = self.settings["newtontol"].GetDouble()  # Tolerance of Newton iterations

        # Initialization
        self.u = np.ones(self.m + 2) * self.ureference  # Velocity
        self.un = np.array(self.u)  # Previous velocity
        self.p = np.ones(self.m + 2) * self.preference  # Kinematic pressure
        self.pn = np.array(self.p)  # Previous kinematic pressure (only value at outlet is used)
        self.a = np.ones(self.m + 2) * np.pi * self.d ** 2 / 4.0  # Area of cross section
        self.an = np.array(self.a)  # Previous area of cross section

        self.disp = np.zeros((self.m, 3))  # Displacement
        self.pres = np.zeros(self.m)  # Pressure
        self.trac = np.zeros((self.m, 3))  # Traction

        self.conditioning = np.pi * self.d ** 2 / 4.0 / (self.ureference + self.dz / self.dt)  # Factor for conditioning Jacobian

        # ModelParts
        self.variable_disp = vars(KM)["DISPLACEMENT"]
        self.variable_pres = vars(KM)["PRESSURE"]
        self.variable_trac = vars(KM)["TRACTION"]
        self.model = cs_data_structure.Model()
        self.model_part = self.model.CreateModelPart("wall")
        self.model_part.AddNodalSolutionStepVariable(self.variable_disp)
        self.model_part.AddNodalSolutionStepVariable(self.variable_pres)
        self.model_part.AddNodalSolutionStepVariable(self.variable_trac)
        for i in range(len(self.z)):
            self.model_part.CreateNewNode(i, 0.0, self.d / 2, self.z[i])
        step = 0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(self.variable_disp, step, self.disp[0, :].tolist())
            node.SetSolutionStepValue(self.variable_pres, step, self.pres[0])
            node.SetSolutionStepValue(self.variable_trac, step, self.trac[0, :].tolist())

        # Interfaces
        self.interface_input = CoSimulationInterface(self.model, self.settings["interface_input"])
        self.interface_output = CoSimulationInterface(self.model, self.settings["interface_output"])

        # Debug
        self.debug = False  # Set on true to save solution of each time step
        self.OutputSolutionStep()

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.k = 0
        self.n += 1
        self.un = np.array(self.u)
        self.pn = np.array(self.p)
        self.an = np.array(self.a)

    def SolveSolutionStep(self, interface_input):
        # Input
        self.disp = interface_input.GetNumpyArray()
        self.interface_input.SetNumpyArray(self.disp)
        a = np.pi * (self.d + 2.0 * self.disp[1::3]) ** 2 / 4.0

        # Input does not contain boundary conditions
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
                if self.inlet_variable == "velocity":
                    self.u[0] = self.GetInletBoundary()
                elif self.inlet_variable == "pressure":
                    self.p[0] = self.GetInletBoundary()
                f = self.GetResidual()
                residual = np.linalg.norm(f)
                if residual / residual0 < self.newtontol:
                    converged = True
                    break
            if not converged:
                Exception("Newton failed to converge")

        self.k += 1
        if self.debug:
            p = self.p[1:self.m + 1] * self.rhof
            u = self.u[1:self.m + 1]
            file_name = self.working_directory + f"/Pressure_Velocity_TS{self.n}_IT{self.k}"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'velocity':<22}\n")
                for i in range(len(self.z)):
                    file.write(f"{self.z[i]:<22}\t{p[i]:<22}\t{u[i]:<22}\n")
            a = self.a[1:self.m + 1]
            file_name = self.working_directory + f"/Area_TS{self.n}_IT{self.k}"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{a[i]:<22}\n')

        # Output does not contain boundary conditions
        self.pres = self.p[1:self.m + 1] * self.rhof
        index_pres = 0
        index_trac = 0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(vars(KM)["PRESSURE"], 0, self.pres[index_pres])
            node.SetSolutionStepValue(vars(KM)["TRACTION"], 0, self.trac[index_trac, :].tolist())
            index_pres += 1
            index_trac += 1
        return self.interface_output.deepcopy()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

    def OutputSolutionStep(self):
        if self.debug:
            p = self.p[1:self.m + 1] * self.rhof
            u = self.u[1:self.m + 1]
            file_name = self.working_directory + f"/Pressure_Velocity_TS{self.n}"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'velocity':<22}\n")
                for i in range(len(self.z)):
                    file.write(f"{self.z[i]:<22}\t{p[i]:<22}\t{u[i]:<22}\n")

    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()

    def SetInterfaceInput(self):
        raise Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()

    def SetInterfaceOutput(self):
        raise Exception("This solver interface provides no mapping.")

    def GetInletBoundary(self):
        if self.inlet_type == 1:
            x = self.inlet_reference \
                + self.inlet_amplitude * np.sin(2.0 * np.pi * (self.n * self.dt) / self.inlet_period)
        elif self.inlet_type == 2:
            x = self.inlet_reference + (self.inlet_amplitude if self.n <= self.inlet_period / self.dt else 0.0)
        elif self.inlet_type == 3:
            x = self.inlet_reference \
                + self.inlet_amplitude * (np.sin(np.pi * (self.n * self.dt) / self.inlet_period)) ** 2
        else:
            x = self.inlet_reference + self.inlet_amplitude * (self.n * self.dt) / self.inlet_period
        return x

    def GetResidual(self):
        usign = self.u[1:self.m + 1] > 0
        ur = self.u[1:self.m + 1] * usign + self.u[2:self.m + 2] * (1.0 - usign)
        ul = self.u[0:self.m] * usign + self.u[1:self.m + 1] * (1.0 - usign)
        self.alpha = np.pi * self.d ** 2 / 4.0 / (self.ureference + self.dz / self.dt)

        f = np.zeros(2 * self.m + 4)
        if self.inlet_variable == "velocity":
            f[0] = (self.u[0] - self.GetInletBoundary()) * self.conditioning
            f[1] = (self.p[0] - (2.0 * self.p[1] - self.p[2])) * self.conditioning
        elif self.inlet_variable == "pressure":
            f[0] = (self.u[0] - (2.0 * self.u[1] - self.u[2])) * self.conditioning
            f[1] = (self.p[0] - self.GetInletBoundary()) * self.conditioning
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
        f[2 * self.m + 2] = (self.u[self.m + 1] - (2.0 * self.u[self.m] - self.u[self.m - 1])) * self.conditioning
        if self.outlet_type == 1:
            f[2 * self.m + 3] = (self.p[self.m + 1] - (2.0 * (self.cmk2 - (np.sqrt(self.cmk2
                                                                                   - self.pn[self.m + 1] / 2.0)
                                                              - (self.u[self.m + 1] - self.un[self.m + 1]) / 4.0) ** 2))
                                 ) * self.conditioning
        else:
            f[2 * self.m + 3] = (self.p[self.m + 1] - self.preference) * self.conditioning
        return f

    def GetJacobian(self):
        usign = self.u[1:self.m + 1] > 0
        j = np.zeros((self.Al + self.Au + 1, 2 * self.m + 4))
        if self.inlet_variable == "velocity":
            j[self.Au + 0 - 0, 0] = 1.0 * self.conditioning  # [0,0]
            j[self.Au + 1 - 1, 1] = 1.0 * self.conditioning  # [1,1]
            j[self.Au + 1 - 3, 3] = -2.0 * self.conditioning  # [1,3]
            j[self.Au + 1 - 5, 5] = 1.0 * self.conditioning  # [1,5]
        elif self.inlet_variable == "pressure":
            j[self.Au + 0 - 0, 0] = 1.0 * self.conditioning  # [0,0]
            j[self.Au + 0 - 2, 2] = -2.0 * self.conditioning  # [0,2]
            j[self.Au + 0 - 4, 4] = 1.0 * self.conditioning  # [0,4]
            j[self.Au + 1 - 1, 1] = 1.0 * self.conditioning  # [1,1]

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
                                              * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0)  # [2*i+1, 2*(i+1)]
        j[self.Au - 3, 5:2 * self.m + 5:2] = -self.alpha  # [2*i, 2*(i+1)+1]
        j[self.Au - 2, 5:2 * self.m + 5:2] = (self.a[1:self.m + 1]
                                              + self.a[2:self.m + 2]) / 4.0  # [2*i+1, 2*(i+1)+1]

        j[self.Au + (2 * self.m + 2) - (2 * self.m + 2), 2 * self.m + 2] = 1.0 * self.conditioning  # [2*m+2, 2*m+2]
        j[self.Au + (2 * self.m + 2) - (2 * self.m), 2 * self.m] = -2.0 * self.conditioning  # [2*m+2, 2*m]
        j[self.Au + (2 * self.m + 2) - (2 * self.m - 2), 2 * self.m - 2] = 1.0 * self.conditioning  # [2*m+2, 2*m-2]
        if self.outlet_type == 1:
            j[self.Au + (2 * self.m + 3) - (2 * self.m + 2), 2 * self.m + 2] = (-(np.sqrt(self.cmk2
                                                                                          - self.pn[self.m + 1] / 2.0)
                                                                                  - (self.u[self.m + 1]
                                                                                     - self.un[self.m + 1])
                                                                                  / 4.0)) * self.conditioning  # [2*m+3, 2*m+2]
        j[self.Au + (2 * self.m + 3) - (2 * self.m + 3), 2 * self.m + 3] = 1.0 * self.conditioning  # [2*m+3, 2*m+3]
        return j
