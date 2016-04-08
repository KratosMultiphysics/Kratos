from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import os
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import DEM_procedures
import shutil

def AddExtraDofs(project_parameters, fluid_model_part, spheres_model_part, cluster_model_part, DEM_inlet_model_part):

    if VELOCITY_LAPLACIAN in project_parameters.fluid_vars:
        for node in fluid_model_part.Nodes:
            node.AddDof(VELOCITY_LAPLACIAN_X)
            node.AddDof(VELOCITY_LAPLACIAN_Y)
            node.AddDof(VELOCITY_LAPLACIAN_Z)

    if VELOCITY_LAPLACIAN_RATE in project_parameters.fluid_vars:
        for node in fluid_model_part.Nodes:
            node.AddDof(VELOCITY_LAPLACIAN_RATE_X)
            node.AddDof(VELOCITY_LAPLACIAN_RATE_Y)
            node.AddDof(VELOCITY_LAPLACIAN_RATE_Z)

def RenumberNodesIdsToAvoidRepeating(fluid_model_part, dem_model_part, rigid_faces_model_part):

    max_id = 1
    renumerate = False

    for node in fluid_model_part.Nodes:

        if node.Id > max_id:
            max_id = node.Id

    for node in dem_model_part.Nodes:

        if node.Id < max_id:
            renumerate = True
            break

    if renumerate:

        print("WARNING!, the DEM model part and the fluid model part have some ID values in common")
        print("Renumerating DEM model part and fem-DEM model parts Ids")

        for node in dem_model_part.Nodes:
            node.Id += max_id

        for node in rigid_faces_model_part.Nodes:
            node.Id += max_id


def RenumberModelPartNodesFromGivenId(model_part, id):

    new_id = id + 1

    for node in model_part.Nodes:
        node.Id = new_id
        new_id = new_id + 1

def RenumberModelPartElementsFromGivenId(model_part, id):

    new_id = id + 1

    for element in model_part.Elements:
        element.Id = new_id
        new_id = new_id + 1


def SetModelPartSolutionStepValue(model_part, var, value):

    for node in model_part.Nodes:
        node.SetSolutionStepValue(var, 0, value)

def InitializeVariablesWithNonZeroValues(fluid_model_part, balls_model_part, pp):
    if pp.CFD_DEM.coupling_level_type:
        SetModelPartSolutionStepValue(fluid_model_part, FLUID_FRACTION, 1.0)
        SetModelPartSolutionStepValue(balls_model_part, FLUID_FRACTION_PROJECTED, 1.0)

def FixModelPart(model_part):

    for node in model_part.Nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_Z)

def GetWordWithSpaces(word, total_length):

    for i in range(len(word), total_length):
        word += ' '

    return word
    
def TransferFacePressuresToPressure(model_part):
    
    for node in model_part.Nodes:
        total_pressure = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE) + node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE)
        node.SetSolutionStepValue(PRESSURE, total_pressure)         

def Norm(my_list):
    return math.sqrt(sum([value ** 2 for value in my_list]))

def FindClosestNode(model_part, coors):
     relative_coors_nodes = [[node.X - coors[0], node.Y - coors[1], node.Z - coors[2]] for node in model_part.Nodes]
     nodes = [node for node in model_part.Nodes]
     min_dist = Norm(coors_nodes[0])
     min_i = 0
     for i in range(len(nodes)):
         norm_i = Norm(relative_coors_nodes[i])
         if min_dist > norm_i:
             min_dist = norm_i
             min_i = i

     return nodes[min_i]

class FluidFractionFieldUtility:

    def CheckIsInside(self, p, l, h):  # p is the point, l and h are the low and high corners of the bounding box

        if p[0] < l[0] or p[1] < l[1] or p[2] < l[2] or p[0] > h[0] or p[1] > h[1] or p[2] > h[2]:
            return False

        else:
            return True

    class LinearField:

        def __init__(self,
                     fluid_fraction_at_zero,
                     fluid_fraction_gradient,
                     lower_corner,
                     upper_corner):

            self.frac_0 = fluid_fraction_at_zero
            self.frac_grad = fluid_fraction_gradient
            self.low = lower_corner
            self.high = upper_corner

    def __init__(self, fluid_model_part, min_fluid_fraction):
        self.fluid_model_part = fluid_model_part
        self.min_fluid_fraction = min_fluid_fraction
        self.field_list = list()

    def AppendLinearField(self, field):
        self.field_list.append(field)

    def AddFluidFractionField(self):

        print('******************************************************************')
        print()
        print('Adding Imposed Fluid Fraction Fields...')
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
                fluid_fraction = node.GetSolutionStepValue(FLUID_FRACTION, 0)

                if self.CheckIsInside([node.X, node.Y, node.Z], field.low, field.high):
                    value = fluid_fraction + field.frac_0 + field.frac_grad[0] * node.X + field.frac_grad[1] * node.Y + field.frac_grad[2] * node.Z
                    value = min(max(value, 0.0), self.min_fluid_fraction)
                    node.SetSolutionStepValue(FLUID_FRACTION, 0, value)

def MultiplyNodalVariableByFactor(model_part, variable, factor):

    for node in model_part.Nodes:
        new_variable = node.GetSolutionStepValue(variable, 0) * factor
        node.SetSolutionStepValue(variable, 0, new_variable)

def ApplySimilarityTransformations(fluid_model_part, transformation_type, mod_over_real):

    if transformation_type == 0:
        return

    elif transformation_type == 1:

        print('***\n\nWARNING!, applying similarity transformations to the problem fluid variables')
        print('The particles diameters quotient is\n')
        print('D_model / D_real =', mod_over_real)
        print()

        if transformation_type == 1:  # Tsuji 2013, (Preserves Archimedes and Reynolds numbers)

            print ('The fluid variables to be modified are\n\nDENSITY\nVISCOSITY\n\n***')

            fluid_density_factor = mod_over_real
            fluid_viscosity_factor = mod_over_real * mod_over_real
            MultiplyNodalVariableByFactor(fluid_model_part, DENSITY, fluid_density_factor)
            MultiplyNodalVariableByFactor(fluid_model_part, VISCOSITY, fluid_viscosity_factor)
    else:

        print(('The entered value similarity_transformation_type = ', transformation_type, 'is not currently supported'))


def FindMaxNodeIdInFLuid(fluid_model_part):

    max = 0

    for node in fluid_model_part.Nodes:

        if node.Id > max:
            max = node.Id

    return max
    
    
def FindMaxElementIdInFLuid(fluid_model_part):

    max = 0

    for element in fluid_model_part.Elements:

        if element.Id > max:
            max = element.Id

    return max    

def FunctionsCalculator(pp):
    if pp.domain_size == 2:
        custom_functions_tool = CustomFunctionsCalculator2D()

    elif pp.domain_size == 3:
        custom_functions_tool = CustomFunctionsCalculator3D()

    return custom_functions_tool

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
        
            shutil.rmtree(main_path + '/' + name, ignore_errors = True)
        
            if not os.path.isdir(dir_abs_path):
                os.makedirs(str(dir_abs_path))

        return directories

    def ControlEcho(self, step, incremental_time, total_steps_expected):

        if incremental_time > self.param.ControlTime:
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

        if estimated_sim_duration / 86400 > 2.0:

            print(('WARNING!!!:       VERY LASTING CALCULATION' + '\n'))

class ProjectionDebugUtils:

    def __init__(self, domain_volume, fluid_model_part, particles_model_part, custom_functions_calculator):
        self.balls_model_part = particles_model_part
        self.fluid_model_part = fluid_model_part
        self.custom_utils = custom_functions_calculator
        self.UpdateDataAndPrint(domain_volume, False)

    def UpdateDataAndPrint(self, domain_volume, is_time_to_print = True):
        self.granul_utils                    = DEM_procedures.GranulometryUtils(domain_volume, self.balls_model_part)
        self.domain_volume                   = domain_volume
        self.number_of_balls                 = self.balls_model_part.NumberOfElements(0)
        self.discr_domain_volume             = self.custom_utils.CalculateDomainVolume(self.fluid_model_part)
        self.proj_fluid_volume               = self.custom_utils.CalculateGlobalFluidVolume(self.fluid_model_part)
        self.solid_volume                    = self.granul_utils.solid_volume
        self.balls_per_area                  = self.granul_utils.spheres_per_area
        self.fluid_volume                    = domain_volume - self.solid_volume
        self.discr_fluid_volume              = self.discr_domain_volume - self.solid_volume
        self.proj_solid_volume               = self.discr_domain_volume - self.proj_fluid_volume
        self.global_fluid_fraction           = self.fluid_volume / self.domain_volume
        self.global_solid_fraction           = 1.0 - self.global_fluid_fraction
        self.fluid_on_balls_total_force      = Array3()
        self.proj_balls_on_fluid_total_force = Array3()
        self.custom_utils.CalculateTotalHydrodynamicForceOnParticles(self.balls_model_part, self.fluid_on_balls_total_force)
        self.custom_utils.CalculateTotalHydrodynamicForceOnFluid(self.fluid_model_part, self.proj_balls_on_fluid_total_force)

        if not is_time_to_print:
            return

        # printing

        tot_len = 32 # total length of each line, including spaces
        print()
        print("Projection-related measurements")
        print("*****************************************************")
        print(GetWordWithSpaces("number_of_balls", tot_len)                 + '=', self.number_of_balls)
        print(GetWordWithSpaces("domain_volume", tot_len)                   + '=', self.domain_volume)
        print(GetWordWithSpaces("fluid_volume", tot_len)                    + '=', self.fluid_volume)
        print(GetWordWithSpaces("solid_volume", tot_len)                    + '=', self.solid_volume)
        print(GetWordWithSpaces("discr_domain_volume", tot_len)             + '=', self.discr_domain_volume)
        print(GetWordWithSpaces("discr_fluid_volume", tot_len)              + '=', self.discr_fluid_volume)
        print(GetWordWithSpaces("proj_fluid_volume", tot_len)               + '=', self.proj_fluid_volume)
        print(GetWordWithSpaces("proj_solid_volume", tot_len)               + '=', self.proj_solid_volume)
        print(GetWordWithSpaces("global_fluid_fraction", tot_len)           + '=', self.global_fluid_fraction)
        print(GetWordWithSpaces("global_solid_fraction", tot_len)           + '=', self.global_solid_fraction)
        print(GetWordWithSpaces("balls_per_area", tot_len)                  + '=', self.balls_per_area)
        print(GetWordWithSpaces("fluid_on_balls_total_force", tot_len)      + '=', self.fluid_on_balls_total_force)
        print(GetWordWithSpaces("proj_balls_on_fluid_total_force", tot_len) + '=', self.proj_balls_on_fluid_total_force)
        print("*****************************************************")
        print()

# This class is useful to keep track of cycles in loops. It is initialized by giving the number of steps per cycle,
# the step at which the cycle starts and weather it is active or not (Tick() returns False in this case).
# Tick() adds 1 to the general counter and to the cycle counter every time it is called, returning True when
# the general counter is greater than the 'beginning_step' and cycling counter is back to the beggining of the cycle
# (and it is active); and False otherwise.

class Counter:

    def __init__(self, steps_in_cycle = 1, beginning_step = 1, is_active = True):

        if steps_in_cycle <= 0 or not isinstance(steps_in_cycle , int):
            raise ValueError("Error: The input steps_in_cycle must be a strictly positive integer")

        self.beginning_step = beginning_step
        self.step           = 1
        self.steps_in_cycle = steps_in_cycle
        self.step_in_cycle  = steps_in_cycle
        self.is_active      = is_active

    def Tick(self):           

        if self.step < self.beginning_step or not self.is_active:
            self.step += 1
            return False

        if self.step_in_cycle == self.steps_in_cycle:
            self.step += 1
            self.step_in_cycle = 1
            return True

        else:
            self.step += 1
            self.step_in_cycle += 1
            return False

    def SetActivation(self, is_active):
        self.is_active = is_active

    def Activate(self, activate):
        self.is_active = self.is_active or activate

    def Deactivate(self, deactivate):
        self.is_active = self.is_active and not deactivate

    def GetStep(self):
        return self.step

    def GetStepInCycle(self):
        return self.step_in_cycle


class PostUtils:

    def __init__(self,
                 gid_io,
                 project_parameters,
                 fluid_model_part,
                 balls_model_part,
                 clusters_model_part,
                 rigid_faces_model_part,
                 mixed_model_part):

        self.gid_io                 = gid_io
        self.fluid_model_part       = fluid_model_part
        self.balls_model_part       = balls_model_part
        self.clusters_model_part    = clusters_model_part
        self.rigid_faces_model_part = rigid_faces_model_part
        self.mixed_model_part       = mixed_model_part
        self.pp                     = project_parameters
        self.post_utilities         = PostUtilities()

    def Writeresults(self, time):

        print("")
        print("*******************  PRINTING RESULTS FOR GID  ***************************")
        sys.stdout.flush()

        if self.pp.GiDMultiFileFlag == "Multiples":
            self.mixed_model_part.Elements.clear()
            self.mixed_model_part.Nodes.clear()
            # here order is important!
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.balls_model_part)
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.rigid_faces_model_part)
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.fluid_model_part)

        self.gid_io.write_swimming_DEM_results(time, self.fluid_model_part, self.balls_model_part, self.clusters_model_part, self.rigid_faces_model_part, self.mixed_model_part, self.pp.nodal_results, self.pp.dem_nodal_results, self.pp.clusters_nodal_results, self.pp.rigid_faces_nodal_results, self.pp.mixed_nodal_results, self.pp.gauss_points_results)

    def ComputeMeanVelocitiesinTrap(self, file_name, time_dem):

        if self.pp.dem.VelocityTrapOption:
            average_velocity = Array3()
            low_point = Array3()
            low_point[0] = self.pp.dem.VelocityTrapMinX
            low_point[1] = self.pp.dem.VelocityTrapMinY
            low_point[2] = self.pp.dem.VelocityTrapMinZ
            high_point = Array3()
            high_point[0] = self.pp.dem.VelocityTrapMaxX
            high_point[1] = self.pp.dem.VelocityTrapMaxY
            high_point[2] = self.pp.dem.VelocityTrapMaxZ

            average_velocity = self.post_utilities.VelocityTrap(self.balls_model_part, low_point, high_point)
            f = open(file_name, 'a')
            tmp = str(time_dem) + "   " + str(average_velocity[0]) + "   " + str(average_velocity[1]) + "   " + str(average_velocity[2]) + "\n"
            f.write(tmp)
            f.flush()
            f.close()

class ResultsFileCreator:
    def __init__(self, model_part, node_id, scalar_vars_list = None, vector_vars_list = None):
        self.file_name = 'results_node_' + str(node_id) + '.txt'

        if scalar_vars_list is None:
            scalar_vars_list = []

        if vector_vars_list is None:
            vector_vars_list = []

        self.scalar_vars = scalar_vars_list

        for var in self.scalar_vars:
            a = str(var).split()
            print(a)

        self.vector_vars = vector_vars_list
        self.n_scalars = len(scalar_vars_list)
        self.n_vectors = len(vector_vars_list)
        self.results = []

    def Record(self, model_part, node_id, time):
        line = [None] * (1 + self.n_scalars + 3 * self.n_vectors)
        line[0] = time

        i = 1
        for var in self.scalar_vars:
            line[i] = model_part.Nodes[node_id].GetSolutionStepValue(var)
            i += 1

        for var in self.vector_vars:
            value = model_part.Nodes[node_id].GetSolutionStepValue(var)
            line[i] = value[0]
            i += 1
            line[i] = value[1]
            i += 1
            line[i] = value[2]
            i += 1

        self.results.append(line)

    def PrintFile(self):
        with open(self.file_name, mode = 'wt') as f:
            header = ['Time']
            for var in self.scalar_vars:
                header.append(' ' + str(var).split()[0])
            for var in self.vector_vars:
                header.append(' ' + str(var).split()[0] + '_X')
                header.append(' ' + str(var).split()[0] + '_Y')
                header.append(' ' + str(var).split()[0] + '_Z')
            for entry in header:
                f.write(entry)
            f.write('\n')

            for result in self.results:
                line = ''
                for entry in result:
                    line += str('%.17f' % entry) + ' '
                f.write(line + ' \n')

class StationarityAssessmentTool:

    def __init__(self, max_pressure_variation_rate_tol, custom_functions_tool):
        self.tol  = max_pressure_variation_rate_tol
        self.tool = custom_functions_tool

    def Assess(self, model_part): # in the first time step the 'old' pressure vector is created and filled
        stationarity = self.tool.AssessStationarity(model_part, self.tol)

        if stationarity:
            print("**************************************************************************************************")
            print()
            print("The model has reached a stationary state. The fluid calculation is suspended.")
            print()
            print("**************************************************************************************************")

        return stationarity
