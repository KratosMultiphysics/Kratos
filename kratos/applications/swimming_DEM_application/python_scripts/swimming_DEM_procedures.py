from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import os
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import DEM_procedures


def AddNodalVariables(model_part, variable_list):

    for var in variable_list:
        model_part.AddNodalSolutionStepVariable(var)

def AddingDEMProcessInfoVariables(pp, dem_model_part):
    dem_model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, pp.buoyancy_force_type)
    dem_model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, pp.drag_force_type)
    dem_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, pp.virtual_mass_force_type)
    dem_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, pp.lift_force_type)
    dem_model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, pp.fluid_model_type)
    dem_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, pp.manually_imposed_drag_law_option)
    dem_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, pp.drag_modifier_type)
    dem_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, pp.initial_drag_force)
    dem_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, pp.drag_law_slope)
    dem_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, pp.power_law_tol)

def ConstructListsOfResultsToPrint(pp):
    pp.dem_nodal_results = []

    if (pp.dem.PostRadius):
        pp.dem_nodal_results += ["RADIUS"]

    if (pp.dem.PostAngularVelocity):
        pp.dem_nodal_results += ["ANGULAR_VELOCITY"]

    if (pp.projection_module_option):

        if (pp.print_REYNOLDS_NUMBER_option):
            pp.dem_nodal_results += ["REYNOLDS_NUMBER"]

        if (pp.print_PRESSURE_GRAD_PROJECTED_option):
            pp.dem_nodal_results += ["PRESSURE_GRAD_PROJECTED"]

        if (pp.print_FLUID_VEL_PROJECTED_option):
            pp.dem_nodal_results += ["FLUID_VEL_PROJECTED"]

    if (pp.print_BUOYANCY_option > 0):
        pp.dem_nodal_results += ["BUOYANCY"]

    if (pp.add_each_hydro_force_option):

        if (pp.drag_force_type > 0 and pp.print_DRAG_FORCE_option):
            pp.dem_nodal_results += ["DRAG_FORCE"]

        if (pp.virtual_mass_force_type > 0 and pp.print_VIRTUAL_MASS_FORCE_option):
            pp.dem_nodal_results += ["VIRTUAL_MASS_FORCE"]

        if (pp.lift_force_type > 0 and pp.print_LIFT_FORCE_option):
            pp.dem_nodal_results += ["LIFT_FORCE"]

    pp.mixed_nodal_results = ["VELOCITY", "DISPLACEMENT"]

    pp.variables_to_print_in_file = ["DRAG_FORCE", "LIFT_FORCE", "BUOYANCY", "VELOCITY"]


def ModifyProjectParameters(pp):
    # defining and adding imposed porosity fields
    pp.solid_fraction_fields = []
    field1 = SolidFractionFieldUtility.LinearField(0.0,
                                                  [0.0, 0.0, - 1.0],
                                                  [-1.0, -1.0, 0.15],
                                                  [1.0, 1.0, 0.3])

    field2 = SolidFractionFieldUtility.LinearField(1.0,
                                                  [0.0, 0.0, 1.0],
                                                  [-1.0, -1.0, 0.3],
                                                  [1.0, 1.0, 0.45])

    pp.solid_fraction_fields.append(field1)
    pp.solid_fraction_fields.append(field2)

    # changes on PROJECT PARAMETERS for the sake of consistency
    ChangeListOfFluidNodalResultsToPrint(pp)

    # applying changes to input data to avoid inconsistencies
    ChangeInputDataForConsistency(pp)

def ChangeListOfFluidNodalResultsToPrint(pp):

    if (pp.coupling_level_type > 0 and pp.print_SOLID_FRACTION_option):
        pp.nodal_results += ["SOLID_FRACTION"]

    if (pp.fluid_model_type == 0 and pp.print_MESH_VELOCITY1_option):
        pp.nodal_results += ["MESH_VELOCITY1"]

    if (pp.fluid_model_type == 1 and pp.print_SOLID_FRACTION_GRADIENT_option):
        pp.nodal_results += ["SOLID_FRACTION_GRADIENT"]

    if (pp.body_force_on_fluid_option and pp.print_BODY_FORCE_option):
        pp.nodal_results += ["BODY_FORCE"]

    if (pp.coupling_level_type > 0 and pp.print_HYDRODYNAMIC_REACTION_option):
        pp.nodal_results += ["HYDRODYNAMIC_REACTION"]

def ChangeInputDataForConsistency(pp):
    pp.dem.project_from_particles_option *= pp.projection_module_option
    pp.project_at_every_substep_option *= pp.projection_module_option

    if (pp.flow_in_porous_medium_option):
        pp.coupling_weighing_type = - 1 # the solid fraction is not projected from DEM (there may not be a DEM part) but externally imposed

    pp.time_steps_per_stationarity_step = max(1, int(pp.time_steps_per_stationarity_step)) # it should never be smaller than 1!
    pp.stationary_problem_option *= not pp.dem.project_from_particles_option

    if (pp.dem.project_from_particles_option or pp.fluid_model_type == 0):
        pp.coupling_level_type = 1

    if (pp.coupling_level_type == 1):
        pp.virtual_mass_force_type = 0

    for var in pp.mixed_nodal_results:

        if var in pp.nodal_results:
            pp.nodal_results.remove(var)

def ConstructListsOfVariables(pp):
    # COUPLING VARIABLES

    # fluid coupling variables
    pp.coupling_fluid_vars = []
    pp.coupling_fluid_vars += [ACCELERATION]

    if (pp.body_force_on_fluid_option):
        pp.coupling_fluid_vars += [BODY_FORCE]

    if (pp.fluid_model_type == 0):
        pp.coupling_fluid_vars += [MESH_VELOCITY1]

    if (pp.fluid_model_type == 0 or pp.coupling_level_type == 1 or pp.drag_force_type == 4):
        pp.coupling_fluid_vars += [SOLID_FRACTION]

    if (pp.fluid_model_type == 1):
        pp.coupling_fluid_vars += [SOLID_FRACTION_GRADIENT]
        pp.coupling_fluid_vars += [SOLID_FRACTION_RATE]

    if (pp.coupling_level_type == 1):
        pp.coupling_fluid_vars += [HYDRODYNAMIC_REACTION]

    # dem coupling variables
    pp.coupling_dem_vars = []

    if (pp.projection_module_option):
        pp.coupling_dem_vars += [FLUID_VEL_PROJECTED]
        pp.coupling_dem_vars += [FLUID_DENSITY_PROJECTED]
        pp.coupling_dem_vars += [PRESSURE_GRAD_PROJECTED]
        pp.coupling_dem_vars += [FLUID_VISCOSITY_PROJECTED]
        pp.coupling_dem_vars += [HYDRODYNAMIC_FORCE]

    if (pp.coupling_level_type == 1 or pp.fluid_model_type == 0):
       pp.coupling_dem_vars += [SOLID_FRACTION_PROJECTED]

    if (pp.lift_force_type == 1):
       pp.coupling_dem_vars += [FLUID_VORTICITY_PROJECTED]
       pp.coupling_dem_vars += [SHEAR_RATE_PROJECTED]

    if (pp.virtual_mass_force_type == 1):
       pp.coupling_dem_vars += [FLUID_ACC_PROJECTED]

    if (pp.embedded_option):
       pp.coupling_dem_vars += [DISTANCE]

    if (pp.drag_force_type == 2):
       pp.coupling_dem_vars += [POWER_LAW_N]
       pp.coupling_dem_vars += [POWER_LAW_K]
       pp.coupling_dem_vars += [GEL_STRENGTH]
       pp.coupling_dem_vars += [YIELD_STRESS]

    # VARIABLES TO ADD
    # listing nodal variables to be added to the model parts (memory will be allocated for them)

    # fluid coupling variables
    pp.fluid_vars = []
    pp.fluid_vars += pp.coupling_fluid_vars
    pp.fluid_vars += [PRESSURE_GRADIENT, AUX_DOUBLE_VAR]

    if (pp.drag_force_type == 2):
        pp.fluid_vars += [POWER_LAW_N]
        pp.fluid_vars += [POWER_LAW_K]
        pp.fluid_vars += [GEL_STRENGTH]
        pp.fluid_vars += [YIELD_STRESS]
        pp.fluid_vars += [BINGHAM_SMOOTHER]

    # dem coupling variables
    pp.dem_vars = []
    pp.dem_vars += pp.coupling_dem_vars

    if (pp.print_REYNOLDS_NUMBER_option):
       pp.dem_vars += [REYNOLDS_NUMBER]

    if (pp.buoyancy_force_type > 0):
       pp.dem_vars += [BUOYANCY]

    if (pp.drag_force_type > 0 and  pp.add_each_hydro_force_option):
       pp.dem_vars += [DRAG_FORCE]

    if (pp.lift_force_type > 0 and  pp.add_each_hydro_force_option):
       pp.dem_vars += [LIFT_FORCE]

    if (pp.virtual_mass_force_type > 0 and  pp.add_each_hydro_force_option):
       pp.dem_vars += [VIRTUAL_FORCE]

    # fem-dem coupling variables
    pp.fem_dem_vars = [VELOCITY, DISPLACEMENT]

    if (pp.embedded_option):
        pp.fem_dem_vars += [FORCE]
        pp.fem_dem_vars += [POSITIVE_FACE_PRESSURE]
        pp.fem_dem_vars += [NEGATIVE_FACE_PRESSURE]

    # inlet coupling variables
    pp.inlet_vars = pp.dem_vars

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

        elif (var == HYDRODYNAMIC_REACTION):

            for node in model_part.Nodes:
                node.SetSolutionStepValue(HYDRODYNAMIC_REACTION_X, 0, 0.0)
                node.SetSolutionStepValue(HYDRODYNAMIC_REACTION_Y, 0, 0.0)
                node.SetSolutionStepValue(HYDRODYNAMIC_REACTION_Z, 0, 0.0)

def FixModelPart(model_part):

    for node in model_part.Nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_Z)

def GetWordWithSpaces(word, total_length):

    for i in range(len(word), total_length):
        word += ' '

    return word

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

    def __init__(self, domain_volume, fluid_model_part, particles_model_part):

        self.balls_model_part = particles_model_part
        self.fluid_model_part = fluid_model_part
        self.UpdateData(domain_volume)

    def UpdateData(self, domain_volume):
        self.granul_utils        = DEM_procedures.GranulometryUtils(domain_volume, self.balls_model_part)
        self.custom_utils        = CustomFunctionsCalculator()
        self.domain_volume       = domain_volume
        self.number_of_balls     = self.balls_model_part.NumberOfElements(0)
        self.discr_domain_volume = self.custom_utils.CalculateDomainVolume(self.fluid_model_part)
        self.proj_fluid_volume   = self.custom_utils.CalculateGlobalFluidVolume(self.fluid_model_part)
        self.solid_volume        = self.granul_utils.solid_volume
        self.fluid_volume        = domain_volume - self.solid_volume
        self.discr_fluid_volume  = self.discr_domain_volume - self.solid_volume
        self.proj_solid_volume   = self.discr_domain_volume - self.proj_fluid_volume
        self.global_porosity     = self.fluid_volume / self.domain_volume
        self.balls_per_area      = domain_volume / self.number_of_balls

    def PrintCurrentData(self):

        tot_len = 25 # total length of the line including spaces
        print()
        print("Volume-related mesurements")
        print("************************************************")
        print(GetWordWithSpaces("number_of_balls", tot_len)     + '=', self.number_of_balls)
        print(GetWordWithSpaces("domain_volume", tot_len)       + '=', self.domain_volume)
        print(GetWordWithSpaces("fluid_volume", tot_len)        + '=', self.fluid_volume)
        print(GetWordWithSpaces("solid_volume", tot_len)        + '=', self.solid_volume)
        print(GetWordWithSpaces("discr_domain_volume", tot_len) + '=', self.discr_domain_volume)
        print(GetWordWithSpaces("discr_fluid_volume", tot_len)  + '=', self.discr_fluid_volume)
        print(GetWordWithSpaces("proj_fluid_volume", tot_len)   + '=', self.proj_fluid_volume)
        print(GetWordWithSpaces("proj_solid_volume", tot_len)   + '=', self.proj_solid_volume)
        print(GetWordWithSpaces("global_porosity", tot_len)     + '=', self.global_porosity)
        print(GetWordWithSpaces("balls_per_area", tot_len)      + '=', self.balls_per_area)
        print("************************************************")
        print()

class ProjectionModule:

    def __init__(self, fluid_model_part, balls_model_part, FEM_DEM_model_part, dimension, pp):

        self.fluid_model_part = fluid_model_part
        self.particles_model_part = balls_model_part
        self.FEM_DEM_model_part = FEM_DEM_model_part
        self.dimension = dimension
        self.max_solid_fraction = pp.max_solid_fraction
        self.coupling_type = pp.coupling_weighing_type
        self.n_particles_in_depth = pp.n_particles_in_depth

        if (self.dimension == 3):
            self.projector = BinBasedDEMFluidCoupledMapping3D(self.max_solid_fraction, self.coupling_type)
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(fluid_model_part)

        else:
            self.projector = BinBasedDEMFluidCoupledMapping2D(self.max_solid_fraction, self.coupling_type, self.n_particles_in_depth)
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(fluid_model_part)

        # telling the projector which variables we are interested in modifying

        for var in pp.coupling_dem_vars:
            self.projector.AddDEMCouplingVariable(var)

        for var in pp.coupling_fluid_vars:
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

class PostUtils:

    def __init__(self,
                 gid_io,
                 project_parameters,
                 fluid_model_part,
                 balls_model_part,
                 fem_dem_model_part,
                 mixed_model_part):

        self.gid_io             = gid_io
        self.fluid_model_part   = fluid_model_part
        self.balls_model_part   = balls_model_part
        self.fem_dem_model_part = fem_dem_model_part
        self.mixed_model_part   = mixed_model_part
        self.pp                 = project_parameters
        self.post_utilities     = PostUtilities()

    def Writeresults(self, time):

        print("")
        print("*******************  PRINTING RESULTS FOR GID  ***************************")
        sys.stdout.flush()

        if (self.pp.GiDMultiFileFlag == "Multiples"):
            self.mixed_model_part.Elements.clear()
            self.mixed_model_part.Nodes.clear()
            # here order is important!
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.balls_model_part)
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.fem_dem_model_part)
            self.post_utilities.AddModelPartToModelPart(self.mixed_model_part, self.fluid_model_part)

        self.gid_io.write_swimming_DEM_results(time, self.fluid_model_part, self.balls_model_part, self.fem_dem_model_part, self.mixed_model_part, self.pp.nodal_results, self.pp.dem_nodal_results, self.pp.mixed_nodal_results, self.pp.gauss_points_results)

    def ComputeMeanVelocitiesinTrap(self, file_name, time_dem):

        if (self.pp.dem.VelocityTrapOption):
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
