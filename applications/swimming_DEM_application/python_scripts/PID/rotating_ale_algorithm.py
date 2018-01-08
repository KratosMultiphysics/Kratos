from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from DEM_procedures import KratosPrint as Say
import pre_calculated_fluid_algorithm
BaseAlgorithm = pre_calculated_fluid_algorithm.Algorithm
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
        self.UU = np.array([a * axis for a in axis])
        self.Ux = np.array([[0, - axis[2], axis[1]],
                            [axis[2], 0., -axis[0]],
                            [-axis[1], axis[0], 0.]])

    def Rotate(self, model_part, time):
        Say('Starting mesh movement...')

        sin = math.sin(self.omega * time)
        cos = math.cos(self.omega * time)

        # Rotation matrix
        R = cos * self.I + sin * self.Ux + (1.0 - cos) * self.UU

        # Rotation matrix' (derivative of R with respect to time)
        Rp = - self.omega * sin * self.I + self.omega * cos * self.Ux + self.omega * sin * self.UU

        for node in model_part.Nodes:
            P0 = np.array([node.X0, node.Y0, node.Z0])

            P = self.a_init + R.dot(P0 - self.a_init)

            Displacement = P - P0
            Velocity = Rp.dot(P0 - self.a_init)

            node.X, node.Y, node.Z = P[0], P[1], P[2]

            node.SetSolutionStepValue(DISPLACEMENT, Vector(list(Displacement)))
            node.SetSolutionStepValue(MESH_VELOCITY, Vector(list(Velocity)))

        Say('Mesh movement finshed.')

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
        self.rotator = MeshRotationUtility(self.pp.CFD_DEM)

    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)
        self.pp.CFD_DEM.AddEmptyValue("ALE_option").SetBool(True)

    def UpdateALEMeshMovement(self, time):
        if self.pp.CFD_DEM["ALE_option"].GetBool():
            self.rotator.RotateMesh(self.fluid_model_part, time)
            self.projection_module.UpdateDatabase(self.h_min)
