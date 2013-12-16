import math
import os
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *

def AddNodalVariables(model_part, variable_list):

    for var in variable_list:
        model_part.AddNodalSolutionStepVariable(var)

def RenumberNodesIdsToAvoidRepeating(fluid_model_part, dem_model_part, fem_dem_model_part):

    max_id = 1
    renumerate = False

    for node in fluid_model_part.Nodes:

        if (node.Id > max_id):
            max_id = node.Id

    for node in dem_model_part.Nodes:

        if (node.Id < max_id):
            renumerate = True
            break

    if (renumerate):
        print "WARNING!, the DEM model part and the fluid model part have some ID values in common"
        print "Renumerating DEM model part and fem-DEM model parts Ids"

        for node in dem_model_part.Nodes:
            node.Id += max_id

        for node in fem_dem_model_part.Nodes:
            node.Id += max_id

def InitializeVariablesToZero(model_part, variable_list):

    for var in variable_list:

        if (var == SOLID_FRACTION):

            for node in model_part.Nodes:
                node.SetSolutionStepValue(var, 0, 0.0)

        elif (var == MESH_VELOCITY1):

            for node in model_part.Nodes:
                node.SetSolutionStepValue(MESH_VELOCITY1_X, 0, 0.0)
                node.SetSolutionStepValue(MESH_VELOCITY1_Y, 0, 0.0)
                node.SetSolutionStepValue(MESH_VELOCITY1_Z, 0, 0.0)

        elif (var == DRAG_REACTION):

            for node in model_part.Nodes:
                node.SetSolutionStepValue(DRAG_REACTION_X, 0, 0.0)
                node.SetSolutionStepValue(DRAG_REACTION_Y, 0, 0.0)
                node.SetSolutionStepValue(DRAG_REACTION_Z, 0, 0.0)

def MultiplyNodalVariableByFactor(model_part, variable, factor):

    for node in model_part.Nodes:
        new_variable = node.GetSolutionStepValue(variable, 0) * factor
        node.SetSolutionStepValue(variable, 0, new_variable)

def ApplySimilarityTransformations(fluid_model_part, transformation_type, mod_over_real):

    if (transformation_type == 0):
        return

    elif (transformation_type == 1):
        print ('***\n\nWARNING!, applying similarity transformations to the problem fluid variables')
        print ('The particles diameters quotient is\n')
        print 'D_model / D_real =', mod_over_real; print

        if (transformation_type == 1): # Tsuji 2013, (Preserves Archimedes and Reynolds numbers)
            print ('The fluid variables to be modified are\n\nDENSITY\nVISCOSITY\n\n***')
            fluid_density_factor = mod_over_real
            fluid_viscosity_factor = mod_over_real * mod_over_real
            MultiplyNodalVariableByFactor(fluid_model_part, DENSITY, fluid_density_factor)
            MultiplyNodalVariableByFactor(fluid_model_part, VISCOSITY, fluid_viscosity_factor)
    else:
        print ('The entered value similarity_transformation_type = ', transformation_type, 'is not currently supported')

def yield_DEM_time(current_time, current_time_plus_increment, delta_time):
    current_time += delta_time

    while current_time < current_time_plus_increment - delta_time:
        yield current_time
        current_time += delta_time

    current_time = current_time_plus_increment
    yield current_time
    

def FindMaxNodeIdInFLuid(fluid_model_part):
    
    max = 0
    for node in fluid_model_part.Nodes:
        if ( node.Id > max ):
            max = node.Id
            
    return max        


class IOTools:

    def __init__(self, Param):

        self.param = Param

    def PrintParticlesResults(self, variables, time, model_part):

        for variablename in variables:
            outstring = "Particles_" + variablename + ".csv"
            outputfile = open(outstring, 'a')
            variables_dictionary = {"PRESSURE" : PRESSURE,
                                    "VELOCITY" : VELOCITY,
                                    "BUOYANCY" : BUOYANCY,
                                    "DRAG_FORCE" : DRAG_FORCE,
                                    "MU" : MU}

            for node in model_part.Nodes:
                Results_value = node.GetSolutionStepValue(variables_dictionary[variablename])
                outputfile.write(str(time) + " " + str(Results_value[0]) + " " +str(Results_value[1]) + " " + str(Results_value[2]) + "\n")

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

        estimated_sim_duration = 60.0 * (total_steps_expected / step) #seconds
        print('The total calculation estimated time is ' + str(estimated_sim_duration) + 'seconds.' + '\n')
        print('In minutes :' + str(estimated_sim_duration / 60) + 'min.'  + '\n')
        print('In hours :' + str(estimated_sim_duration / 3600) + 'hrs.'  + '\n')
        print('In days :' + str(estimated_sim_duration / 86400) + 'days.' + '\n')

        if (estimated_sim_duration / 86400 > 2.0):
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

    def __init__(self, fluid_model_part, balls_model_part, FEM_DEM_model_part, dimension, max_solid_fraction, coupling_type, n_particles_in_depth):

        self.fluid_model_part     = fluid_model_part
        self.particles_model_part = balls_model_part
        self.FEM_DEM_model_part   = FEM_DEM_model_part
        self.dimension            = dimension
        self.max_solid_fraction   = max_solid_fraction
        self.coupling_type        = coupling_type
        self.n_particles_in_depth = n_particles_in_depth

        if (self.dimension == 3):
            self.projector = BinBasedDEMFluidCoupledMapping3D(max_solid_fraction, coupling_type)
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(fluid_model_part)

        else:
            self.projector = BinBasedDEMFluidCoupledMapping2D(max_solid_fraction, coupling_type, n_particles_in_depth)
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(fluid_model_part)

        # calculating the fluid nodal areas that are needed for the coupling

        self.area_calculator = CalculateNodalAreaProcess(self.fluid_model_part, self.dimension)
        self.area_calculator.Execute()

    def UpdateDatabase(self, HMin):

        if (self.dimension == 3):
            self.bin_of_objects_fluid.UpdateSearchDatabase()

        else:
            self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ProjectFromFluid(self, alpha):

        self.projector.InterpolateFromFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid, alpha)

    def ProjectFromNewestFluid(self):

        self.projector.InterpolateFromNewestFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ProjectFromParticles(self):

        self.projector.InterpolateFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.bin_of_objects_fluid)

    def ComputePostProcessResults(self, particles_process_info):
        self.projector.ComputePostProcessResults(self.particles_model_part, self.fluid_model_part, self.FEM_DEM_model_part, self.bin_of_objects_fluid, particles_process_info)

