from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import os
from KratosMultiphysics import *
#from KratosMultiphysics.IncompressibleFluidApplication import *
#from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import DEM_procedures
import shutil
import os
import weakref

def Say(*args):
    Logger.PrintInfo("DEM-FLUID", *args)
    Logger.Flush()

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

    max_fluid_id = FindMaxNodeId(fluid_model_part)
    must_renumber = True in (node.Id < max_fluid_id for node in dem_model_part.Nodes + rigid_faces_model_part.Nodes)

    if must_renumber:

        Logger.PrintWarning("DEM-FLUID","WARNING!, the DEM model part and the fluid model part have some ID values in common. Renumbering...")

        for node in dem_model_part.Nodes:
            node.Id += max_fluid_id

        for node in rigid_faces_model_part.Nodes:
            node.Id += max_fluid_id

        Logger.PrintWarning("DEM-FLUID","The DEM model part and the fem-DEM model parts Ids have been renumbered")


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
    checker = VariableChecker()

    if checker.ModelPartHasNodalVariableOrNot(fluid_model_part, FLUID_FRACTION):
        SetModelPartSolutionStepValue(fluid_model_part, FLUID_FRACTION, 1.0)
    if checker.ModelPartHasNodalVariableOrNot(balls_model_part, FLUID_FRACTION_PROJECTED):
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

def NormOfDifference(v1, v2):
    return math.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)

def Norm(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

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

        Logger.PrintInfo("DEM-FLUID",'******************************************************************')
        Logger.PrintInfo()
        Logger.PrintInfo("DEM-FLUID",'Adding Imposed Fluid Fraction Fields...')
        Logger.PrintInfo()
        Logger.Flush()

        count = 0

        for field in self.field_list:
            count += 1

            Logger.PrintInfo("DEM-FLUID",'field number', count, ':')
            Logger.PrintInfo()
            Logger.PrintInfo("DEM-FLUID",vars(field))
            Logger.PrintInfo()

        Logger.PrintInfo("DEM-FLUID",'******************************************************************')
        Logger.Flush()

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

        Logger.PrintWarning("DEM-FLUID",'***\n\nWARNING!, applying similarity transformations to the problem fluid variables')
        Logger.PrintWarning("DEM-FLUID",'The particles diameters quotient is\n')
        Logger.PrintWarning("DEM-FLUID",'D_model / D_real =', mod_over_real)
        Logger.PrintWarning("DEM-FLUID",)

        if transformation_type == 1:  # Tsuji 2013, (Preserves Archimedes and Reynolds numbers)

            Logger.PrintWarning ('The fluid variables to be modified are\n\nDENSITY\nVISCOSITY\n\n***')

            fluid_density_factor = mod_over_real
            fluid_viscosity_factor = mod_over_real * mod_over_real
            MultiplyNodalVariableByFactor(fluid_model_part, DENSITY, fluid_density_factor)
            MultiplyNodalVariableByFactor(fluid_model_part, VISCOSITY, fluid_viscosity_factor)
    else:

        Logger.PrintWarning("DEM-FLUID",('The entered value similarity_transformation_type = ', transformation_type, 'is not currently supported'))


def FindMaxNodeId(fluid_model_part):
    return max((node.Id for node in fluid_model_part.Nodes))

def FindMaxElementId(fluid_model_part):
    return max((element.Id for element in fluid_model_part.Elements))

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

            Say('Real time calculation: ' + str(incremental_time))
            Say('Percentage Completed: ' + str(percentage) + ' %')
            Say("TIME STEP = " + str(step) + '\n')

            prev_time = (incremental_time)

class ProjectionDebugUtils:

    def __init__(self, domain_volume, fluid_model_part, particles_model_part, custom_functions_calculator):
        self.balls_model_part = particles_model_part
        self.fluid_model_part = fluid_model_part
        self.custom_utils = custom_functions_calculator
        self.UpdateDataAndPrint(domain_volume, False)

    def UpdateDataAndPrint(self, domain_volume, is_time_to_print = True):
        self.granul_utils                         = DEM_procedures.GranulometryUtils(domain_volume, self.balls_model_part)
        self.domain_volume                        = domain_volume
        self.number_of_balls                      = self.balls_model_part.NumberOfElements(0)
        self.discr_domain_volume                  = self.custom_utils.CalculateDomainVolume(self.fluid_model_part)
        self.proj_fluid_volume                    = self.custom_utils.CalculateGlobalFluidVolume(self.fluid_model_part)
        self.solid_volume                         = self.granul_utils.solid_volume
        self.balls_per_area                       = self.granul_utils.spheres_per_area
        self.fluid_volume                         = domain_volume - self.solid_volume
        self.discr_fluid_volume                   = self.discr_domain_volume - self.solid_volume
        self.proj_solid_volume                    = self.discr_domain_volume - self.proj_fluid_volume
        self.global_fluid_fraction                = self.fluid_volume / self.domain_volume
        self.global_solid_fraction                = 1.0 - self.global_fluid_fraction
        self.fluid_on_balls_total_force           = Array3()
        self.proj_balls_on_fluid_total_force      = Array3()
        self.mean_proj_balls_on_fluid_total_force = Array3()
        self.custom_utils.CalculateTotalHydrodynamicForceOnParticles(self.balls_model_part, self.fluid_on_balls_total_force)
        self.custom_utils.CalculateTotalHydrodynamicForceOnFluid(self.fluid_model_part, self.proj_balls_on_fluid_total_force, self.mean_proj_balls_on_fluid_total_force)

        if not is_time_to_print:
            return

        # printing

        tot_len = 38 # total length of each line, including spaces
        Logger.PrintInfo("DEM-FLUID",)
        Logger.PrintInfo("DEM-FLUID","Projection-related measurements")
        Logger.PrintInfo("DEM-FLUID",tot_len * "**")
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("number_of_balls", tot_len)                      + '=', self.number_of_balls)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("domain_volume", tot_len)                        + '=', self.domain_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("fluid_volume", tot_len)                         + '=', self.fluid_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("solid_volume", tot_len)                         + '=', self.solid_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("discr_domain_volume", tot_len)                  + '=', self.discr_domain_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("discr_fluid_volume", tot_len)                   + '=', self.discr_fluid_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("proj_fluid_volume", tot_len)                    + '=', self.proj_fluid_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("proj_solid_volume", tot_len)                    + '=', self.proj_solid_volume)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("global_fluid_fraction", tot_len)                + '=', self.global_fluid_fraction)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("global_solid_fraction", tot_len)                + '=', self.global_solid_fraction)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("balls_per_area", tot_len)                       + '=', self.balls_per_area)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("fluid_on_balls_total_force", tot_len)           + '=', self.fluid_on_balls_total_force)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("proj_balls_on_fluid_total_force", tot_len)      + '=', self.proj_balls_on_fluid_total_force)
        Logger.PrintInfo("DEM-FLUID",GetWordWithSpaces("mean_proj_balls_on_fluid_total_force", tot_len) + '=', self.mean_proj_balls_on_fluid_total_force)
        Logger.PrintInfo("DEM-FLUID",tot_len * "**")
        Logger.PrintInfo("DEM-FLUID",)
        Logger.Flush()

# This class is useful to keep track of cycles in loops. It is initialized by giving the number of steps per cycle,
# the step at which the cycle starts and weather it is active or not (Tick() returns False in this case).
# Tick() adds 1 to the general counter and to the cycle counter every time it is called, returning True when
# the general counter is greater than the 'beginning_step' and cycling counter is back to the beggining of the cycle
# (and it is active); and False otherwise.

class Counter:

    def __init__(self,
                 steps_in_cycle = 1,
                 beginning_step = 1,
                 is_active = True,
                 is_dead = False):

        if steps_in_cycle <= 0 or not isinstance(steps_in_cycle , int):
            raise ValueError("Error: The input steps_in_cycle must be a strictly positive integer")

        self.beginning_step    = beginning_step
        self.step              = 1
        self.steps_in_cycle    = steps_in_cycle
        self.step_in_cycle     = steps_in_cycle
        self.is_active         = is_active
        self.is_dead           = is_dead
        self.accumulated_ticks = 0

    def Tick(self):

        if self.is_dead:
            return False

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

    def SuperTick(self, n_ticks = 1):
        self.accumulated_ticks += 1
        if self.accumulated_ticks >= n_ticks:
            self.accumulated_ticks = 0
            return True
        else:
            return False

    def SetActivation(self, is_active):
        self.is_active = is_active

    def Activate(self, condition = True):
        self.is_active |= condition

    def Deactivate(self, condition = True):
        self.is_active &= not condition

    def Switch(self, condition = None):
        if condition == None:
            self.is_active = not self.is_active
        else:
            self.is_active = condition

    def GetStep(self):
        return self.step

    def GetStepInCycle(self):
        return self.step_in_cycle

    def Kill(self):
        self.is_dead = True




class Averager:
    def __init__(self, steps_in_cycle = 1, beginning_step = 1, is_active = True):
        self.counter = Counter(steps_in_cycle, beginning_step, is_active)
        self.min = float('+inf')
        self.max = float('-inf')
        self.sum = 0.0
        self.sum_error = 0.0
        self.average = 0.0
        self.average_error = 0.0
        self.step = 0
    def Norm(self, v):
        if self.counter.Tick():
            self.step += 1
            value = math.sqrt(sum([entry ** 2 for entry in v]))
            self.sum += value
            self.min = min(value, self.min)
            self.max = max(value, self.max)
            self.average = self.sum / self.step
    def RelativeError(self, v_exact, v):
        if self.counter.Tick():
            self.step += 1
            norm_exact = math.sqrt(sum([entry ** 2 for entry in v_exact]))
            self.sum_exact += norm_exact
            self.average_exact = self.sum / self.step
            value = math.sqrt(sum([(v_exact[i] - v[i]) ** 2 / self.average_exact for entry in v]))
            self.sum += value
            self.min = min(value, self.min)
            self.max = max(value, self.max)
            self.average_error = self.sum_error / self.step

    def GetCurrentData(self):
        return self.min, self.max, self.average

class PostUtils:

    def __init__(self,
                 gid_io,
                 project_parameters,
                 fluid_model_part,
                 balls_model_part,
                 clusters_model_part,
                 rigid_faces_model_part,
                 mixed_model_part):

        self.gid_io                 = weakref.proxy(gid_io)
        self.fluid_model_part       = fluid_model_part
        self.balls_model_part       = balls_model_part
        self.clusters_model_part    = clusters_model_part
        self.rigid_faces_model_part = rigid_faces_model_part
        self.mixed_model_part       = mixed_model_part
        self.pp                     = project_parameters
        self.post_utilities         = PostUtilities()

    def Writeresults(self, time):

        Logger.PrintInfo("DEM-FLUID","")
        Logger.PrintInfo("DEM-FLUID","*******************  PRINTING RESULTS FOR GID  ***************************")
        Logger.Flush()

        if self.pp.GiDMultiFileFlag == "Multiples":
            self.mixed_model_part.Elements.clear()
            self.mixed_model_part.Nodes.clear()
            # here order is important!
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.balls_model_part)
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.rigid_faces_model_part)
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.fluid_model_part)

        self.gid_io.write_swimming_DEM_results(time,
                                               self.fluid_model_part,
                                               self.balls_model_part,
                                               self.clusters_model_part,
                                               self.rigid_faces_model_part,
                                               self.mixed_model_part,
                                               self.pp.nodal_results,
                                               self.pp.dem_nodal_results,
                                               self.pp.clusters_nodal_results,
                                               self.pp.rigid_faces_nodal_results,
                                               self.pp.mixed_nodal_results,
                                               self.pp.gauss_points_results)

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

# The following function creates a run_code to be appended to the name of the PostFiles directory for the benchmark marine_rain (2013 Guseva)
def CreateRunCode(pp):
    code = []

    if pp.CFD_DEM["basset_force_type"].GetInt() > 0:
        history_or_not = 'H'
    else:
        history_or_not = 'NH'

    code.append(history_or_not)

    if pp.CFD_DEM["basset_force_type"].GetInt() == 4:
        method_name = 'Hinsberg'
        number_of_exponentials = 'm=' + str(pp.CFD_DEM.number_of_exponentials)
        time_window = 'tw=' + str(pp.CFD_DEM["time_window"].GetDouble())
        code.append(method_name)
        code.append(number_of_exponentials)
        code.append(time_window)

    elif pp.CFD_DEM["basset_force_type"].GetInt() > 0:
        method_name = 'Daitche'
        code.append(method_name)
    else:
        method_name = pp.CFD_DEM["TranslationalIntegrationScheme"].GetString()
        code.append(method_name)

    DEM_dt = 'Dt=' + str(pp.CFD_DEM["MaxTimeStep"].GetDouble())
    code.append(DEM_dt)

    if pp.CFD_DEM["basset_force_type"].GetInt() > 0:
        phi = 'phi=' + str(round(1 / pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt(), 3))
        code.append(phi)

    quadrature_order = 'QuadOrder=' + str(pp.CFD_DEM["quadrature_order"].GetInt())
    code.append(quadrature_order)

    return '_' + '_'.join(code)

def CopyInputFilesIntoFolder(files_path, folder_path):
    import glob, os, shutil
    input_data_directory = folder_path + "/InputData"
    if not os.path.isdir(input_data_directory):
        os.makedirs(str(input_data_directory))

    file_endings = ["*.mdpa", "DEM_explicit_solver_var.py", "ProjectParameters.py", "marine_rain.py", "*.geo", "*.dat", "*.msh"]
    for ending in file_endings:
        files = glob.iglob(os.path.join(files_path, ending))

        for file in files:
            if os.path.isfile(file):
                shutil.copy2(file, input_data_directory)

def PrintPositionsToFile(coors, method='', j_var='', k_var='', main_path=''):

    with open(main_path + '/coors_' + str(method) + str(j_var) + str(k_var) + '.txt','w') as fout:
        for coor in coors:
            line = ''

            for comp in coor:
                line += str(comp) + ' '
            fout.write(line + '\n')

def PrintPositionsToFileInterpolation(coors, n_div, method, main_path=''):
    if method == 2:
        mat_deriv_type = 'recovery'
    else:
        mat_deriv_type = 'standard'

    directory_name = main_path + '/coordinates'
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    with open(main_path + '/coordinates/coors_ndiv_' + str(n_div) + '_type_' + mat_deriv_type + '.txt','w') as fout:
        for coor in coors:
            line = ''

            for comp in coor:
                line += str(comp) + ' '
            fout.write(line + '\n')

class StationarityAssessmentTool:

    def __init__(self, max_pressure_variation_rate_tol, custom_functions_tool):
        self.tol  = max_pressure_variation_rate_tol
        self.tool = weakref.proxy(custom_functions_tool)

    def Assess(self, model_part): # in the first time step the 'old' pressure vector is created and filled
        stationarity = self.tool.AssessStationarity(model_part, self.tol)

        if stationarity:
            Logger.PrintInfo("DEM-FLUID","**************************************************************************************************")
            Logger.PrintInfo()
            Logger.PrintInfo("DEM-FLUID","The model has reached a stationary state. The fluid calculation is suspended.")
            Logger.PrintInfo()
            Logger.PrintInfo("DEM-FLUID","**************************************************************************************************")
            Logger.Flush()

        return stationarity
