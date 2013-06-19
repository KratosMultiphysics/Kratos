import math

from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

def CalculateGlobalPorosity(ParticleModelPart, DomainVolume):
    n_balls      = 0
    solid_volume = 0.0

    for ball in ParticleModelPart.Nodes:
        n_balls += 1
        radius = ball.GetSolutionStepValue(RADIUS)
        volume = 4 / 3 * pow(radius, 3) * math.pi
        solid_volume += volume
        global_porosity = solid_volume / (solid_volume + DomainVolume)

    balls_per_area = DomainVolume / n_balls
    output = [solid_volume, global_porosity, n_balls, balls_per_area]
    return output
