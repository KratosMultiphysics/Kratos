#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  Fluid Dynamics Application
#  License: BSD License (kratos/license.txt)
#  Main author: Nicolò Antonelli
#
import math

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication as CFD
from KratosMultiphysics.FluidDynamicsApplication.cylinder_sbm_geometry import CylinderSbmGeometry
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_lifecycle import MovingGapLifecycle
from KratosMultiphysics.FluidDynamicsApplication.patch_history_transfer import transfer_background_histories


class CylinderSbmGeometryTest(KratosUnittest.TestCase):
    def test_projection_and_translation(self):
        model = KM.Model()
        part = model.CreateModelPart('Patch')
        part.AddNodalSolutionStepVariable(KM.VELOCITY)
        part.SetBufferSize(3)
        properties = part.CreateNewProperties(1)
        properties[KM.CONSTITUTIVE_LAW] = CFD.Newtonian2DLaw()
        for j in range(7):
            for i in range(7):
                part.CreateNewNode(1+j*7+i, i, j, 0.0)
        for j in range(6):
            for i in range(6):
                n = 1+j*7+i
                element = part.CreateNewElement('Element2D4N', 1+j*6+i,
                                                [n, n+1, n+8, n+7], properties)
                element.Set(KM.ACTIVE, True)
        cylinder = CylinderSbmGeometry(part, (3.0, 3.0), 1.3, 0.1, 0.0)
        self.assertGreater(sum(e.IsNot(KM.ACTIVE) for e in part.Elements), 0)
        conditions = part.GetSubModelPart('CylinderBoundary').Conditions
        self.assertGreater(len(conditions), 0)
        cylinder.update((.25, -.1), (.5, -.2))
        for condition in conditions:
            projection = condition.GetValue(KM.SURROGATE_BOUNDARY_PROJECTION)
            face = condition.GetValue(KM.SURROGATE_BOUNDARY_FACE_COORDINATES)
            velocity = condition.GetValue(CFD.SBM_BOUNDARY_VELOCITIES)
            for i in range(2):
                self.assertAlmostEqual(math.hypot(projection[i, 0]-3.25, projection[i, 1]-2.9), 1.3)
                self.assertAlmostEqual(velocity[i, 0], .5)
                self.assertAlmostEqual(velocity[i, 1], -.2)
            dx, dy = face[1, 0]-face[0, 0], face[1, 1]-face[0, 1]
            midpoint = ((face[0, 0]+face[1, 0])/2, (face[0, 1]+face[1, 1])/2)
            self.assertGreater(dy*(3.25-midpoint[0])-dx*(2.9-midpoint[1]), 0.0)
        target_part = model.CreateModelPart('Targets')
        target_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        target_part.SetBufferSize(3)
        target = target_part.CreateNewNode(1000, 3.0, 3.0, 0.0)
        with self.assertRaisesRegex(RuntimeError, 'Missing 1 required velocity histories'):
            transfer_background_histories(part, [target], .1, .1, (0., 0.), required_steps=(1,))

    def test_sinusoidal_historical_displacement(self):
        lifecycle = MovingGapLifecycle.__new__(MovingGapLifecycle)
        lifecycle._initial_time = 0.0
        lifecycle.motions = [{'type': 'sinusoidal', 'amplitude': .1, 'frequency': 1.0}]
        lifecycle.velocities = [(0.0, 0.0)]
        self.assertAlmostEqual(lifecycle._translation(0, .25)[0], .1)
        self.assertAlmostEqual(lifecycle._translation(0, .75)[0], -.1)


if __name__ == '__main__':
    KratosUnittest.main()
