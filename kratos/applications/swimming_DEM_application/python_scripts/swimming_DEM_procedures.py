import math

from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

class PorosityUtils:

    def __init__(self, DomainVolume, ParticlesModelPart):

        self.balls_model_part = ParticlesModelPart
        self.UpdateData(DomainVolume)

    def UpdateData(self, DomainVolume):

        self.number_of_balls = 0
        self.solid_volume    = 0.0

        for ball in self.balls_model_part.Nodes:
            self.number_of_balls += 1
            radius = ball.GetSolutionStepValue(RADIUS)
            volume = 4 / 3 * pow(radius, 3) * math.pi
            self.solid_volume += volume
            self.global_porosity = self.solid_volume / (self.solid_volume + DomainVolume)

        self.balls_per_area = DomainVolume / self.number_of_balls
        self.fluid_volume   = DomainVolume - self.solid_volume

    def PrintCurrentData(self):

        print "solid volume: ", self.solid_volume
        print "global porosity: ", self.global_porosity
        print "number_of_balls: ", self.number_of_balls
        print "balls per area unit: ", self.balls_per_area

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
