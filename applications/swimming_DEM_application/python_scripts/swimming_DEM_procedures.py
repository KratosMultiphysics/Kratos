import math
import os
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *

class IOTools:

    def __init__(self, Param):

        self.param = Param

    def CreateProblemDirectories(self, main_path, dir_names):

        directories = []
        n_dir = len(dir_names)

        for i in range(n_dir):
            directories.append(str(main_path) + '/' + str(self.param.DEM_problem_name) + '_' + dir_names[i])

        for directory in directories:
            if (not os.path.isdir(directory)):
                os.makedirs(str(directory))

        return directories

    def ControlEcho(self, step, incremental_time, total_steps_expected):

        if (incremental_time > self.param.ControlTime):
            percentage = 100.0 * (float(step) / total_steps_expected)
            print 'Real time calculation: ' + str(incremental_time)
            print 'Percentage Completed: '  + str(percentage) + ' %'
            print "TIME STEP = "            + str(step)   + '\n'
            prev_time = (incremental_time)

    def CalculationLengthEstimationEcho(self, step, incremental_time, total_steps_expected):

        estimated_time = 60.0 * (total_steps_expected / step) #seconds
        print('The total calculation estimated time is ' + str(estimation_time) + 'seconds.' + '\n')
        print('In minutes :' + str(estimation_time / 60) + 'min.'  + '\n')
        print('In hours :' + str(estimation_time / 3600) + 'hrs.'  + '\n')
        print('In days :' + str(estimation_time / 86400) + 'days.' + '\n')

        if (estimated_time / 86400 > 2.0):
            print('WARNING!!!:       VERY LASTING CALCULATION'+'\n')

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

    def __init__(self, FluidModelPart, ParticlesModelPart, Dimension, NParticlesInDepth):

        self.fluid_model_part     = FluidModelPart
        self.particles_model_part = ParticlesModelPart
        self.n_particles_in_depth = NParticlesInDepth
        self.dimension            = Dimension

        if (self.dimension == 3):
            self.projector = BinBasedDEMFluidCoupledMapping3D()
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(FluidModelPart)

        else:
            self.projector = BinBasedDEMFluidCoupledMapping2D()
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(FluidModelPart)

    def UpdateDatabase(self, HMin):

        if (self.dimension == 3):
            self.bin_of_objects_fluid.UpdateSearchDatabase()

        else:

            self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ProjectFromFluid(self):

        self.projector.InterpolationFromFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ProjectFromParticles(self):

        self.projector.InterpolationFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.n_particles_in_depth, self.bin_of_objects_fluid)

