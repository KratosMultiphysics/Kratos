from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import t_junction_algorithm
BaseAlgorithm = t_junction_algorithm.Algorithm
import h5py
import numpy as np
import math

class Rotator:
    def __init__(self, 
                 rotation_axis_initial_point,
                 rotation_axis_final_point,
                 angular_velocity_module):
        self.a_init = Vector(rotation_axis_initial_point)
        self.a_final = Vector(rotation_axis_final_point)
        self.omega = angular_velocity_module
        self.axis = Vector(self.Normalize(self.a_final - self.a_init))
        self.CalculateRodriguesMatrices(self.axis)

    def CalculateRodriguesMatrices(self, axis):
        self.I = np.identity(3)
        self.UU = np.array([[a * axis] for a in axis])
        self.Ux = np.array([[0, - axis[2], axis[1]], 
                             [axis[2], 0., -axis[0]], 
                             [-axis[1], axis[0], 0.]])

    def Rotate(self, node, time):
        sin = math.sin(self.omega * time)
        cos = math.cos(self.omega * time)

        # Rotation matrix
        R = cos * self.I + sin * self.Ux + (1.0 - cos) * self.UU

        # Rotation matrix' (derivative of R with respect to time)
        Rp = - self.omega * sin * self.I + self.omega * cos * self.Ux + self.omega * sin * self.UU

        P0 = Vector([node.X0, node.Y0, node.Z0])
        P = Vector(list(self.a_init + R.dot(P0 - self.a_init)))

        Displacement = Vector(P - P0)
        Velocity = Vector(list(Rp.dot(P0 - self.a_init)))

        node.X, node.Y, node.Z = P[0], P[1], P[2]

        node.SetSolutionStepValue(MESH_VELOCITY, Velocity)

        node.SetSolutionStepValue(DISPLACEMENT, Displacement)

    def Normalize(self, v):
        mod_2 = sum([x ** 2 for x in v])

        if mod_2 == 0:
            return v
        else:
            mod_inv = 1.0 / math.sqrt(mod_2)
            return mod_inv * v

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)
        self.SetRotator()

    def SetRotator(self):
        a_init = self.pp.CFD_DEM.AddEmptyValue("frame_rotation_axis_initial_point").GetVector()
        a_final = self.pp.CFD_DEM.AddEmptyValue("frame_rotation_axis_final_point").GetVector()
        omega = self.pp.CFD_DEM.AddEmptyValue("angular_velocity_magnitude").GetDouble()
        self.rotator = Rotator(a_init, a_final, omega)

    def UpdateMeshMovement(self, time):
      
        for node in self.fluid_model_part.Nodes:
            self.rotator.Rotate(node, time)

    def FluidSolve(self, time = 'None', solve_system = True):
        self.UpdateMeshMovement(time)
        BaseAlgorithm.FluidSolve(self, time, solve_system)
    
    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)
        self.pp.CFD_DEM.AddEmptyValue("ALE_option").SetBool(True)
