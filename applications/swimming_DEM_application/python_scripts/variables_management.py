from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import DEM_explicit_solver_var as DEM_parameters

def AddNodalVariables(model_part, variable_list):

    for var in variable_list:
        model_part.AddNodalSolutionStepVariable(var)

def AddingDEMProcessInfoVariables(pp, dem_model_part):

    dem_model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, pp.swim.buoyancy_force_type)
    dem_model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, pp.swim.drag_force_type)
    dem_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, pp.swim.virtual_mass_force_type)
    dem_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, pp.swim.lift_force_type)
    dem_model_part.ProcessInfo.SetValue(MAGNUS_FORCE_TYPE, pp.swim.magnus_force_type)
    dem_model_part.ProcessInfo.SetValue(HYDRO_TORQUE_TYPE, pp.swim.hydro_torque_type)
    dem_model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, pp.swim.fluid_model_type)
    dem_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, pp.swim.manually_imposed_drag_law_option)
    dem_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, pp.swim.drag_modifier_type)
    dem_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, pp.swim.initial_drag_force)
    dem_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, pp.swim.drag_law_slope)
    dem_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, pp.swim.power_law_tol)
    dem_model_part.ProcessInfo.SetValue(DRAG_POROSITY_CORRECTION_TYPE, pp.swim.drag_porosity_correction_type)

# constructing lists of variables to add
# * Performing modifications to the input parameters for concistency (provisional until interface does it)
# * Choosing the variables to be printed
# * Choosing the variables to be passed as a parameter to the constructor of a ProjectionModule
#       instance to be filled with the other phase's info through the coupling process
# * Listing nodal variables to be added to the model parts (memory will be allocated for them).
#       Note that additional variables may be added as well by the fluid and/or DEM strategies.

def ConstructListsOfVariables(pp):

    # INPUT CHANGES FOR CONSISTENCY
    # performing some extra changes on Project Parameters, ensuring consistency in the input
    ChangeInputDataForConsistency(pp)

    # PRINTING VARIABLES
    # constructing lists of variables to be printed
    ConstructListsOfResultsToPrint(pp)

    # COUPLING VARIABLES
    # listing the variables involved in the fluid-particles coupling
    ConstructListsOfVariablesForCoupling(pp)

    # VARIABLES TO ADD
    # listing nodal variables to be added to the model parts (memory will be allocated for them)

    # fluid variables
    pp.fluid_vars = []
    pp.fluid_vars += pp.fluid_printing_vars
    pp.fluid_vars += pp.coupling_fluid_vars
    pp.fluid_vars += [PRESSURE_GRADIENT]

    if (pp.swim.drag_force_type >= 0):
        pp.fluid_vars += [POWER_LAW_N]
        pp.fluid_vars += [POWER_LAW_K]
        #pp.fluid_vars += [GEL_STRENGTH]
        pp.fluid_vars += [YIELD_STRESS]
        pp.fluid_vars += [BINGHAM_SMOOTHER]

    # dem variables
    pp.dem_vars = []
    pp.dem_vars += pp.dem_printing_vars
    pp.dem_vars += pp.coupling_dem_vars

    if (pp.swim.buoyancy_force_type > 0):
        pp.dem_vars += [BUOYANCY]

    if (pp.swim.drag_force_type > 0 and  pp.swim.add_each_hydro_force_option):
        pp.dem_vars += [DRAG_FORCE]

    if (pp.swim.drag_force_type > 1):
        pp.dem_vars += [PARTICLE_SPHERICITY]

    if (pp.swim.lift_force_type > 0 and  pp.swim.add_each_hydro_force_option):
        pp.dem_vars += [LIFT_FORCE]

    if (pp.swim.virtual_mass_force_type > 0 and  pp.swim.add_each_hydro_force_option):
        pp.dem_vars += [VIRTUAL_MASS_FORCE]

    # clusters variables
    pp.clusters_vars = []

    # rigid faces variables
    pp.rigid_faces_vars = [VELOCITY, DISPLACEMENT, ELASTIC_FORCES, PRESSURE, TANGENTIAL_ELASTIC_FORCES, SHEAR_STRESS, NODAL_AREA]

    if (pp.swim.embedded_option):
        pp.rigid_faces_vars += [FORCE]
        pp.rigid_faces_vars += [POSITIVE_FACE_PRESSURE]
        pp.rigid_faces_vars += [NEGATIVE_FACE_PRESSURE]

    # inlet variables
    pp.inlet_vars = pp.dem_vars
   
def ConstructListsOfResultsToPrint(pp):
    pp.dem_nodal_results = []
    pp.clusters_nodal_results = []
    pp.rigid_faces_nodal_results = []

    if (DEM_parameters.PostRadius):
        pp.dem_nodal_results += ["RADIUS"]

    if (DEM_parameters.PostAngularVelocity):
        pp.dem_nodal_results += ["ANGULAR_VELOCITY"]

    if (pp.swim.projection_module_option):

        if (pp.swim.print_REYNOLDS_NUMBER_option):
            pp.dem_nodal_results += ["REYNOLDS_NUMBER"]

        if (pp.swim.print_PRESSURE_GRAD_PROJECTED_option):
            pp.dem_nodal_results += ["PRESSURE_GRAD_PROJECTED"]

        if (pp.swim.print_HYDRODYNAMIC_FORCE_option):
            pp.dem_nodal_results += ["HYDRODYNAMIC_FORCE"]

        if (pp.swim.print_HYDRODYNAMIC_MOMENT_option):
            pp.dem_nodal_results += ["HYDRODYNAMIC_MOMENT"]

        if (pp.swim.print_FLUID_VEL_PROJECTED_option):
            pp.dem_nodal_results += ["FLUID_VEL_PROJECTED"]

        if (pp.swim.print_FLUID_ACCEL_PROJECTED_option):
            pp.dem_nodal_results += ["FLUID_ACCEL_PROJECTED"]            

        if (pp.swim.print_FLUID_FRACTION_PROJECTED_option):
            pp.dem_nodal_results += ["FLUID_FRACTION_PROJECTED"]

        if (pp.swim.print_FLUID_VISCOSITY_PROJECTED_option):
            pp.dem_nodal_results += ["FLUID_VISCOSITY_PROJECTED"]

        if (pp.swim.print_BUOYANCY_option > 0):
            pp.dem_nodal_results += ["BUOYANCY"]

    if (pp.swim.add_each_hydro_force_option):

        if (pp.swim.print_DRAG_FORCE_option):
            pp.dem_nodal_results += ["DRAG_FORCE"]

        if (pp.swim.print_VIRTUAL_MASS_FORCE_option):
            pp.dem_nodal_results += ["VIRTUAL_MASS_FORCE"]

        if (pp.swim.print_LIFT_FORCE_option):
            pp.dem_nodal_results += ["LIFT_FORCE"]
    
    pp.rigid_faces_nodal_results += ["POSITIVE_FACE_PRESSURE"]
    pp.rigid_faces_nodal_results += ["NEGATIVE_FACE_PRESSURE"]

    # changes on the fluid variables to print for the sake of consistency
    ChangeListOfFluidNodalResultsToPrint(pp)
    
    pp.mixed_nodal_results = ["VELOCITY", "DISPLACEMENT"]

    pp.variables_to_print_in_file = ["DRAG_FORCE", "LIFT_FORCE", "BUOYANCY", "VELOCITY"]

    pp.dem_printing_vars = []
    pp.clusters_printing_vars = []
    pp.fluid_printing_vars = []
    pp.rigid_faces_printing_vars = []

    for variable in pp.nodal_results:
        pp.fluid_printing_vars += [eval(variable)]

    for variable in pp.dem_nodal_results:
        pp.dem_printing_vars += [eval(variable)]

    for variable in pp.clusters_nodal_results:
        pp.clusters_printing_vars += [eval(variable)]
        
    for variable in pp.rigid_faces_nodal_results:
        pp.rigid_faces_printing_vars += [eval(variable)]

    for variable in pp.mixed_nodal_results:
        pp.dem_printing_vars += [eval(variable)]
        pp.fluid_printing_vars += [eval(variable)]
        
    for var in pp.mixed_nodal_results:

        if var in pp.nodal_results:
            pp.nodal_results.remove(var)
            
    if (not pp.swim.print_PRESSURE_option):
        pp.nodal_results.remove("PRESSURE")

    EliminateRepeatedValuesFromList(pp.nodal_results)
    EliminateRepeatedValuesFromList(pp.dem_nodal_results)
    EliminateRepeatedValuesFromList(pp.mixed_nodal_results)

def ConstructListsOfVariablesForCoupling(pp):

    # fluid coupling variables
    pp.coupling_fluid_vars = []
    pp.coupling_fluid_vars += [ACCELERATION]

    if (pp.swim.embedded_option):
       pp.coupling_fluid_vars += [DISTANCE]

    pp.coupling_fluid_vars += [BODY_FORCE]

    if (pp.swim.fluid_model_type == 0):
        pp.coupling_fluid_vars += [AVERAGED_FLUID_VELOCITY]

    if (pp.swim.fluid_model_type == 0 or pp.swim.coupling_level_type == 1 or pp.swim.drag_force_type == 4):
        pp.coupling_fluid_vars += [FLUID_FRACTION]
        
        if (pp.swim.print_SOLID_FRACTION_option):
            pp.coupling_fluid_vars += [SOLID_FRACTION]

    if (pp.swim.fluid_model_type == 1):
        pp.coupling_fluid_vars += [FLUID_FRACTION_GRADIENT]
        pp.coupling_fluid_vars += [FLUID_FRACTION_RATE]

    if (pp.swim.coupling_level_type == 1):
        pp.coupling_fluid_vars += [HYDRODYNAMIC_REACTION]

    if (pp.swim.drag_force_type >= 2):
        pp.coupling_fluid_vars += [POWER_LAW_N]
        pp.coupling_fluid_vars += [POWER_LAW_K]
        #pp.coupling_fluid_vars += [GEL_STRENGTH]
        pp.coupling_fluid_vars += [YIELD_STRESS]

    # dem coupling variables
    pp.coupling_dem_vars = []

    if (pp.swim.projection_module_option):
        pp.coupling_dem_vars += [FLUID_VEL_PROJECTED]
        pp.coupling_dem_vars += [FLUID_DENSITY_PROJECTED]
        pp.coupling_dem_vars += [PRESSURE_GRAD_PROJECTED]
        pp.coupling_dem_vars += [FLUID_VISCOSITY_PROJECTED]
        pp.coupling_dem_vars += [HYDRODYNAMIC_FORCE]
        pp.coupling_dem_vars += [HYDRODYNAMIC_MOMENT]

    if (pp.swim.coupling_level_type == 1 or pp.swim.fluid_model_type == 0):
        pp.coupling_dem_vars += [FLUID_FRACTION_PROJECTED]

    if (pp.swim.lift_force_type == 1):
        pp.coupling_dem_vars += [FLUID_VORTICITY_PROJECTED]
        pp.coupling_dem_vars += [SHEAR_RATE_PROJECTED]

    if (pp.swim.virtual_mass_force_type >= 1):
        pp.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]

    if (pp.swim.embedded_option):
        pp.coupling_dem_vars += [DISTANCE]

    if (pp.swim.drag_force_type >= 2):
        pp.coupling_dem_vars += [POWER_LAW_N]
        pp.coupling_dem_vars += [POWER_LAW_K]
        #pp.coupling_dem_vars += [GEL_STRENGTH]
        pp.coupling_dem_vars += [YIELD_STRESS]

    if (pp.swim.print_REYNOLDS_NUMBER_option):
        pp.coupling_dem_vars += [REYNOLDS_NUMBER]

def ChangeListOfFluidNodalResultsToPrint(pp):

    if (pp.swim.print_FLUID_FRACTION_option):
        pp.nodal_results += ["FLUID_FRACTION"]

    if (pp.swim.print_SOLID_FRACTION_option):
        pp.nodal_results += ["SOLID_FRACTION"]

    if (pp.swim.fluid_model_type == 0 and pp.swim.print_AVERAGED_FLUID_VELOCITY_option):
        pp.nodal_results += ["AVERAGED_FLUID_VELOCITY"]

    if (pp.swim.fluid_model_type == 1 and pp.swim.print_FLUID_FRACTION_GRADIENT_option):
        pp.nodal_results += ["FLUID_FRACTION_GRADIENT"]

    if (pp.swim.body_force_on_fluid_option and pp.swim.print_BODY_FORCE_option):
        pp.nodal_results += ["BODY_FORCE"]

    if (pp.swim.print_HYDRODYNAMIC_REACTION_option):
        pp.nodal_results += ["HYDRODYNAMIC_REACTION"]

    if (pp.swim.embedded_option):
        pp.nodal_results += ["DISTANCE"]

def ChangeInputDataForConsistency(pp):
    pp.swim.coupling_level_type = DEM_parameters.project_from_particles_option
    DEM_parameters.project_from_particles_option *= pp.swim.projection_module_option
    pp.swim.project_at_every_substep_option *= pp.swim.projection_module_option
    pp.swim.lift_force_type *= DEM_parameters.consider_lift_force_option
    pp.swim.drag_modifier_type = DEM_parameters.drag_modifier_type
    
    if (DEM_parameters.VelocityTrapOption == "ON"):
        pp.velocity_trap_option = 1
    else: 
        pp.velocity_trap_option = 0

    if (pp.swim.flow_in_porous_medium_option):
        pp.coupling_weighing_type = - 1 # the fluid fraction is not projected from DEM (there may not be a DEM part) but is externally imposed

    pp.swim.time_steps_per_stationarity_step = max(1, int(pp.swim.time_steps_per_stationarity_step)) # it should never be smaller than 1!
    pp.swim.stationary_problem_option *= not DEM_parameters.project_from_particles_option

def EliminateRepeatedValuesFromList(redundant_list):
    clean_list = []

    for var in redundant_list:

        if var in clean_list:
            redundant_list.remove(var)

        clean_list += [var]
