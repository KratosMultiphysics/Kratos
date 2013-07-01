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

class ProjectionModule:

    def __init__(self, FluidModelPart, ParticlesModelPart, DomainSize, NParticlesInDepth):

        self.fluid_model_part     = FluidModelPart
        self.particles_model_part = ParticlesModelPart
        self.n_particles_in_depth = NParticlesInDepth

        if (DomainSize == 3):
            self.projector = BinBasedDEMFluidCoupledMapping3D()
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(FluidModelPart)

        else:
            self.projector = BinBasedDEMFluidCoupledMapping2D()
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(FluidModelPart)

    def UpdateDatabase(self, HMin):
        self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ProjectFromFluid(self):
        self.projector.InterpolationFromFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ProjectFromParticles(self):
        self.projector.InterpolationFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.n_particles_in_depth, self.bin_of_objects_fluid)
