from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *


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

# constructing lists of variables to add
# * Performing modifications to the input parameters for concistency (provisional until interface does it)
# * Choosing the variables to be printed
# * Choosing the variables to be passed as a parameter to the constructor of a ProjectionModule
#       instance to be filled with the other phase's info through the coupling process
# * Listing nodal variables to be added to the model parts (memory will be allocated for them).
#       Note that additional variables may be added as well by the fluid and/or DEM strategies.

def ConstructListsOfVariables(pp):

    # PRINTING VARIABLES
    # constructing lists of variables to be printed
    ConstructListsOfResultsToPrint(pp)

    # INPUT CHANGES FOR CONCISTENCY
    # performing some extra changes on Project Parameters, ensuring consistency in the input
    ModifyProjectParameters(pp)

    # COUPLING VARIABLES
    # listing the variables involved in the fluid-particles coupling
    ConstructListsOfVariablesForCoupling(pp)

    # VARIABLES TO ADD
    # listing nodal variables to be added to the model parts (memory will be allocated for them)

    # fluid variables
    pp.fluid_vars = []
    pp.fluid_vars += pp.fluid_printing_vars
    pp.fluid_vars += pp.coupling_fluid_vars
    pp.fluid_vars += [PRESSURE_GRADIENT, AUX_DOUBLE_VAR]

    if (pp.drag_force_type >= 0):
        pp.fluid_vars += [POWER_LAW_N]
        pp.fluid_vars += [POWER_LAW_K]
        pp.fluid_vars += [GEL_STRENGTH]
        pp.fluid_vars += [YIELD_STRESS]
        pp.fluid_vars += [BINGHAM_SMOOTHER]

    # dem variables
    pp.dem_vars = []
    pp.dem_vars += [SOLID_FRACTION]
    pp.dem_vars += pp.dem_printing_vars
    pp.dem_vars += pp.coupling_dem_vars

    if (pp.buoyancy_force_type > 0):
       pp.dem_vars += [BUOYANCY]

    if (pp.drag_force_type > 0 and  pp.add_each_hydro_force_option):
       pp.dem_vars += [DRAG_FORCE]

    if (pp.drag_force_type == 2 or pp.drag_force_type == 3 and pp.add_each_hydro_force_option):
       pp.dem_vars += [PARTICLE_SPHERICITY]

    if (pp.lift_force_type > 0 and  pp.add_each_hydro_force_option):
       pp.dem_vars += [LIFT_FORCE]

    if (pp.virtual_mass_force_type > 0 and  pp.add_each_hydro_force_option):
       pp.dem_vars += [VIRTUAL_FORCE]

    # fem-dem variables
    pp.fem_dem_vars = [VELOCITY, DISPLACEMENT]

    if (pp.embedded_option):
        pp.fem_dem_vars += [FORCE]
        pp.fem_dem_vars += [POSITIVE_FACE_PRESSURE]
        pp.fem_dem_vars += [NEGATIVE_FACE_PRESSURE]

    # inlet variables
    pp.inlet_vars = pp.dem_vars

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

        if (pp.print_HYDRODYNAMIC_FORCE_option):
            pp.dem_nodal_results += ["HYDRODYNAMIC_FORCE"]

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

    pp.dem_printing_vars = []
    pp.fluid_printing_vars = []

    for variable in pp.nodal_results:
        pp.fluid_printing_vars += [eval(variable)]

    for variable in pp.dem_nodal_results:
        pp.dem_printing_vars += [eval(variable)]

    for variable in pp.mixed_nodal_results:
        pp.dem_printing_vars += [eval(variable)]
        pp.fluid_printing_vars += [eval(variable)]

def ConstructListsOfVariablesForCoupling(pp):

    # fluid coupling variables
    pp.coupling_fluid_vars = []
    pp.coupling_fluid_vars += [ACCELERATION]

    if (pp.embedded_option):
       pp.coupling_fluid_vars += [DISTANCE]

    pp.coupling_fluid_vars += [BODY_FORCE]

    if (pp.fluid_model_type == 0):
        pp.coupling_fluid_vars += [MESH_VELOCITY1]

    if (pp.fluid_model_type == 0 or pp.coupling_level_type == 1 or pp.drag_force_type == 4):
        pp.coupling_fluid_vars += [FLUID_FRACTION]

    if (pp.fluid_model_type == 1):
        pp.coupling_fluid_vars += [FLUID_FRACTION_GRADIENT]
        pp.coupling_fluid_vars += [FLUID_FRACTION_RATE]

    if (pp.coupling_level_type == 1):
        pp.coupling_fluid_vars += [HYDRODYNAMIC_REACTION]

    if (pp.drag_force_type >= 2):
        pp.coupling_fluid_vars += [POWER_LAW_N]
        pp.coupling_fluid_vars += [POWER_LAW_K]
        pp.coupling_fluid_vars += [GEL_STRENGTH]
        pp.coupling_fluid_vars += [YIELD_STRESS]

    # dem coupling variables
    pp.coupling_dem_vars = []

    if (pp.projection_module_option):
        pp.coupling_dem_vars += [FLUID_VEL_PROJECTED]
        pp.coupling_dem_vars += [FLUID_DENSITY_PROJECTED]
        pp.coupling_dem_vars += [PRESSURE_GRAD_PROJECTED]
        pp.coupling_dem_vars += [FLUID_VISCOSITY_PROJECTED]
        pp.coupling_dem_vars += [HYDRODYNAMIC_FORCE]

    if (pp.coupling_level_type == 1 or pp.fluid_model_type == 0):
       pp.coupling_dem_vars += [FLUID_FRACTION_PROJECTED]

    if (pp.lift_force_type == 1):
       pp.coupling_dem_vars += [FLUID_VORTICITY_PROJECTED]
       pp.coupling_dem_vars += [SHEAR_RATE_PROJECTED]

    if (pp.virtual_mass_force_type == 1):
       pp.coupling_dem_vars += [FLUID_ACC_PROJECTED]

    if (pp.embedded_option):
       pp.coupling_dem_vars += [DISTANCE]

    if (pp.drag_force_type >= 2):
       pp.coupling_dem_vars += [POWER_LAW_N]
       pp.coupling_dem_vars += [POWER_LAW_K]
       pp.coupling_dem_vars += [GEL_STRENGTH]
       pp.coupling_dem_vars += [YIELD_STRESS]

def ModifyProjectParameters(pp):

    # changes on PROJECT PARAMETERS for the sake of consistency
    ChangeListOfFluidNodalResultsToPrint(pp)

    # applying changes to input data to avoid inconsistencies
    ChangeInputDataForConsistency(pp)

def ChangeListOfFluidNodalResultsToPrint(pp):

    if (pp.coupling_level_type > 0 and pp.print_FLUID_FRACTION_option):
        pp.nodal_results += ["FLUID_FRACTION"]

    if (pp.print_SOLID_FRACTION_option):
        pp.nodal_results += ["SOLID_FRACTION"]

    if (pp.fluid_model_type == 0 and pp.print_MESH_VELOCITY1_option):
        pp.nodal_results += ["MESH_VELOCITY1"]

    if (pp.fluid_model_type == 1 and pp.print_FLUID_FRACTION_GRADIENT_option):
        pp.nodal_results += ["FLUID_FRACTION_GRADIENT"]

    if (pp.body_force_on_fluid_option and pp.print_BODY_FORCE_option):
        pp.nodal_results += ["BODY_FORCE"]

    if (pp.coupling_level_type > 0 and pp.print_HYDRODYNAMIC_REACTION_option):
        pp.nodal_results += ["HYDRODYNAMIC_REACTION"]

    if (pp.embedded_option):
        pp.nodal_results += ["DISTANCE"]

def ChangeInputDataForConsistency(pp):
    pp.dem.project_from_particles_option *= pp.projection_module_option
    pp.project_at_every_substep_option *= pp.projection_module_option

    if (pp.flow_in_porous_medium_option):
        pp.coupling_weighing_type = - 1 # the fluid fraction is not projected from DEM (there may not be a DEM part) but is externally imposed

    pp.time_steps_per_stationarity_step = max(1, int(pp.time_steps_per_stationarity_step)) # it should never be smaller than 1!
    pp.stationary_problem_option *= not pp.dem.project_from_particles_option

    if (pp.dem.project_from_particles_option or pp.fluid_model_type == 0):
        pp.coupling_level_type = 1

    if (pp.coupling_level_type == 1):
        pp.virtual_mass_force_type = 0

    #if (pp.coupling_level_type > 0 and not pp.body_force_on_fluid_option):
        #pp.coupling_level_type = 0
        #print('! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')
        #print()
        #print('Two-way coupling cannot be performed if body_force_on_fluid_option = 0 !')
        #print('Running one-way copupling simulation...')
        #print()
        #print('! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ')

    for var in pp.mixed_nodal_results:

        if var in pp.nodal_results:
            pp.nodal_results.remove(var)
