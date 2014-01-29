from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
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

        print("WARNING!, the DEM model part and the fluid model part have some ID values in common")
        print("Renumerating DEM model part and fem-DEM model parts Ids")

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


def FixModelPart(model_part):

    for node in model_part.Nodes:
        node.Fix(VELOCITY)


class SolidFractionFieldUtility:

    def CheckIsInside(self, p, l, h):  # p is the point, l and h are the low and high corners of the bounding box

        if (p[0] < l[0] or p[1] < l[1] or p[2] < l[2] or p[0] > h[0] or p[1] > h[1] or p[2] > h[2]):
            return False

        else:
            return True

    class LinearField:

        def __init__(self,
                     solid_fraction_at_zero,
                     solid_fraction_gradient,
                     lower_corner,
                     upper_corner):

            self.frac_0 = solid_fraction_at_zero
            self.frac_grad = solid_fraction_gradient
            self.low = lower_corner
            self.high = upper_corner

    def __init__(self, fluid_model_part, max_solid_fraction):
        self.fluid_model_part = fluid_model_part
        self.max_solid_fraction = max_solid_fraction
        self.field_list = list()

    def AppendLinearField(self, field):
        self.field_list.append(field)

    def AddSolidFractionField(self):

        print('******************************************************************')
        print()
        print('Adding Imposed Solid Fraction Fields...')
        print()

        count = 0

        for field in self.field_list:
            count += 1

            print('field number', count, ':')
            print()
            print(vars(field))
            print()

        print('******************************************************************')

        for field in self.field_list:

            for node in self.fluid_model_part.Nodes:
                solid_fraction = node.GetSolutionStepValue(SOLID_FRACTION, 0)

                if (self.CheckIsInside([node.X, node.Y, node.Z], field.low, field.high)):
                    value = solid_fraction + field.frac_0 + field.frac_grad[0] * node.X + field.frac_grad[1] * node.Y + field.frac_grad[2] * node.Z
                    value = min(max(value, 0.0), self.max_solid_fraction)
                    node.SetSolutionStepValue(SOLID_FRACTION, 0, value)


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
        print('D_model / D_real =', mod_over_real)
        print()

        if (transformation_type == 1):  # Tsuji 2013, (Preserves Archimedes and Reynolds numbers)

            print ('The fluid variables to be modified are\n\nDENSITY\nVISCOSITY\n\n***')

            fluid_density_factor = mod_over_real
            fluid_viscosity_factor = mod_over_real * mod_over_real
            MultiplyNodalVariableByFactor(fluid_model_part, DENSITY, fluid_density_factor)
            MultiplyNodalVariableByFactor(fluid_model_part, VISCOSITY, fluid_viscosity_factor)
    else:

        print(('The entered value similarity_transformation_type = ', transformation_type, 'is not currently supported'))


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

        if (node.Id > max):
            max = node.Id

    return max


class IOTools:

    def __init__(self, Param):

        self.param = Param

    def PrintParticlesResults(self, variables, time, model_part):

        for variablename in variables:
            outstring = "Particles_" + variablename + ".csv"
            outputfile = open(outstring, 'a')
            variables_dictionary = {"PRESSURE": PRESSURE,
                                    "VELOCITY": VELOCITY,
                                    "BUOYANCY": BUOYANCY,
                                    "DRAG_FORCE": DRAG_FORCE,
                                    "MU": MU}

            for node in model_part.Nodes:
                Results_value = node.GetSolutionStepValue(variables_dictionary[variablename])
                outputfile.write(str(time) + " " + str(Results_value[0]) + " " + str(Results_value[1]) + " " + str(Results_value[2]) + "\n")

    def CreateProblemDirectories(self, main_path, dir_names):

        directories = {}

        for name in dir_names:
            dir_abs_path = main_path + '/' + name
            directories[name] = dir_abs_path

            if (not os.path.isdir(dir_abs_path)):
                os.makedirs(str(dir_abs_path))

        return directories

    def ControlEcho(self, step, incremental_time, total_steps_expected):

        if (incremental_time > self.param.ControlTime):
            percentage = 100.0 * (float(step) / total_steps_expected)

            print('Real time calculation: ' + str(incremental_time))
            print('Percentage Completed: ' + str(percentage) + ' %')
            print("TIME STEP = " + str(step) + '\n')

            prev_time = (incremental_time)

    def CalculationLengthEstimationEcho(self, step, incremental_time, total_steps_expected):

        estimated_sim_duration = 60.0 * (total_steps_expected / step)  # seconds

        print(('The total calculation estimated time is ' + str(estimated_sim_duration) + 'seconds.' + '\n'))
        print(('In minutes :' + str(estimated_sim_duration / 60) + 'min.' + '\n'))
        print(('In hours :' + str(estimated_sim_duration / 3600) + 'hrs.' + '\n'))
        print(('In days :' + str(estimated_sim_duration / 86400) + 'days.' + '\n'))

        if (estimated_sim_duration / 86400 > 2.0):

            print(('WARNING!!!:       VERY LASTING CALCULATION' + '\n'))


class PorosityUtils:

    def __init__(self, DomainVolume, ParticlesModelPart):

        self.balls_model_part = ParticlesModelPart
        self.UpdateData(DomainVolume)

    def UpdateData(self, DomainVolume):

        self.number_of_balls = 0
        self.solid_volume = 0.0

        for ball in self.balls_model_part.Nodes:
            self.number_of_balls += 1
            radius = ball.GetSolutionStepValue(RADIUS)
            volume = 4 / 3 * pow(radius, 3) * math.pi
            self.solid_volume += volume
            self.global_porosity = self.solid_volume / (self.solid_volume + DomainVolume)

        self.balls_per_area = DomainVolume / self.number_of_balls
        self.fluid_volume = DomainVolume - self.solid_volume

    def PrintCurrentData(self):

        print("solid volume: ", self.solid_volume)
        print("global porosity: ", self.global_porosity)
        print("number_of_balls: ", self.number_of_balls)
        print("balls per area unit: ", self.balls_per_area)


class ProjectionModule:

    def __init__(self, fluid_model_part, balls_model_part, FEM_DEM_model_part, dimension, max_solid_fraction, coupling_type, n_particles_in_depth, dem_vars, fluid_vars):

        self.fluid_model_part = fluid_model_part
        self.particles_model_part = balls_model_part
        self.FEM_DEM_model_part = FEM_DEM_model_part
        self.dimension = dimension
        self.max_solid_fraction = max_solid_fraction
        self.coupling_type = coupling_type
        self.n_particles_in_depth = n_particles_in_depth

        if (self.dimension == 3):
            self.projector = BinBasedDEMFluidCoupledMapping3D(max_solid_fraction, coupling_type)
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(fluid_model_part)

        else:
            self.projector = BinBasedDEMFluidCoupledMapping2D(max_solid_fraction, coupling_type, n_particles_in_depth)
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(fluid_model_part)

        # telling the projector which variables we are interested in modifying

        for var in dem_vars:
            self.projector.AddDEMCouplingVariable(var)

        for var in fluid_vars:
            self.projector.AddFluidCouplingVariable(var)

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
