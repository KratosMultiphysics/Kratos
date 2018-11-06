from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
#from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

def AddNodalVariables(model_part, variable_list):

    for var in variable_list:
        model_part.AddNodalSolutionStepVariable(var)

def AddingExtraProcessInfoVariables(pp, fluid_model_part, dem_model_part): #DEPRECATED!

    AddExtraProcessInfoVariablesToFluidModelPart(pp, fluid_model_part)
    AddExtraProcessInfoVariablesToDispersePhaseModelPart(pp, dem_model_part)

# constructing lists of variables to add
# * Performing modifications to the input parameters for consistency (provisional until interface does it)
# * Choosing the variables to be printed
# * Choosing the variables to be passed as a parameter to the constructor of a ProjectionModule
#       instance to be filled with the other phase's info through the coupling process
# * Listing nodal variables to be added to the model parts (memory will be allocated for them).
#       Note that additional variables may be added as well by the fluid and/or DEM strategies.

def AddFrameOfReferenceRelatedVariables(pp, model_part):
    frame_of_reference_type = pp.CFD_DEM["frame_of_reference_type"].GetInt()
    model_part.ProcessInfo.SetValue(FRAME_OF_REFERENCE_TYPE, frame_of_reference_type)

    if frame_of_reference_type == 1: # Rotating frame
        angular_velocity_of_frame = Vector(3)
        angular_velocity_of_frame[:] = [pp.CFD_DEM["angular_velocity_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]

        model_part.ProcessInfo.SetValue(ANGULAR_VELOCITY_MOVING_FRAME, angular_velocity_of_frame)

        if frame_of_reference_type >= 2: # Gemeral frame
            angular_velocity_of_frame_old = Vector(3)
            angular_velocity_of_frame_old[:] = [pp.CFD_DEM["angular_velocity_of_frame_old" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
            acceleration_of_frame_origin = Vector(3)
            acceleration_of_frame_origin[:] = [pp.CFD_DEM["acceleration_of_frame_origin" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
            angular_acceleration_of_frame = Vector(3)
            angular_acceleration_of_frame[:] = [pp.CFD_DEM["angular_acceleration_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
            model_part.ProcessInfo.SetValue(ANGULAR_VELOCITY_MOVING_FRAME_OLD, angular_velocity_of_frame_old)
            model_part.ProcessInfo.SetValue(ACCELERATION_MOVING_FRAME_ORIGIN, acceleration_of_frame_origin)
            model_part.ProcessInfo.SetValue(ANGULAR_ACCELERATION_MOVING_FRAME, angular_acceleration_of_frame)

def AddExtraProcessInfoVariablesToFluidModelPart(pp, fluid_model_part):

    AddFrameOfReferenceRelatedVariables(pp, fluid_model_part)

    fluid_model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1)
    gravity = Vector(3)
    if pp.CFD_DEM["body_force_on_fluid_option"].GetBool():
        gravity[0] = pp.CFD_DEM["GravityX"].GetDouble()
        gravity[1] = pp.CFD_DEM["GravityY"].GetDouble()
        gravity[2] = pp.CFD_DEM["GravityZ"].GetDouble()
    fluid_model_part.ProcessInfo.SetValue(GRAVITY, gravity)

    if pp.CFD_DEM["laplacian_calculation_type"].GetInt() == 3: # recovery through solving a system
        fluid_model_part.ProcessInfo.SetValue(COMPUTE_LUMPED_MASS_MATRIX, 1)
    elif pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 4 or pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 5 or pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 6: # recovery through solving a system
        fluid_model_part.ProcessInfo.SetValue(COMPUTE_LUMPED_MASS_MATRIX, 0)

    if pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 5 or pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 6:
         fluid_model_part.ProcessInfo.SetValue(CURRENT_COMPONENT, 0)

    if pp.CFD_DEM["non_newtonian_option"].GetBool():
        fluid_model_part.ProcessInfo.SetValue(YIELD_STRESS, pp.CFD_DEM["yield_stress"].GetDouble())
        fluid_model_part.ProcessInfo.SetValue(REGULARIZATION_COEFFICIENT, pp.CFD_DEM["regularization_coefficient"].GetDouble())
        fluid_model_part.ProcessInfo.SetValue(POWER_LAW_K, pp.CFD_DEM["power_law_k"].GetDouble())
        fluid_model_part.ProcessInfo.SetValue(POWER_LAW_N, pp.CFD_DEM["power_law_n"].GetDouble())

def AddExtraProcessInfoVariablesToDispersePhaseModelPart(pp, dem_model_part):

    AddFrameOfReferenceRelatedVariables(pp, dem_model_part)

    dem_model_part.ProcessInfo.SetValue(COUPLING_TYPE, pp.CFD_DEM["coupling_level_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, pp.CFD_DEM["buoyancy_force_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, pp.CFD_DEM["drag_force_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, pp.CFD_DEM["virtual_mass_force_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(BASSET_FORCE_TYPE, pp.CFD_DEM["basset_force_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, pp.CFD_DEM["lift_force_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(MAGNUS_FORCE_TYPE, pp.CFD_DEM["magnus_force_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(HYDRO_TORQUE_TYPE, pp.CFD_DEM["hydro_torque_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, pp.CFD_DEM["fluid_model_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, pp.CFD_DEM["manually_imposed_drag_law_option"].GetBool())
    dem_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, pp.CFD_DEM["drag_modifier_type"].GetInt())
    dem_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, pp.CFD_DEM["initial_drag_force"].GetDouble())
    dem_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, pp.CFD_DEM["drag_law_slope"].GetDouble())
    dem_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, pp.CFD_DEM["power_law_tol"].GetDouble())
    dem_model_part.ProcessInfo.SetValue(DRAG_POROSITY_CORRECTION_TYPE, pp.CFD_DEM["drag_porosity_correction_type"].GetInt())

    if pp.CFD_DEM["basset_force_type"].GetInt() > 0:
        dem_model_part.ProcessInfo.SetValue(NUMBER_OF_INIT_BASSET_STEPS, pp.CFD_DEM["n_init_basset_steps"].GetInt())
        dem_model_part.ProcessInfo.SetValue(TIME_STEPS_PER_QUADRATURE_STEP, pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt())
        dem_model_part.ProcessInfo.SetValue(LAST_TIME_APPENDING, 0.0)
        dem_model_part.ProcessInfo.SetValue(QUADRATURE_ORDER, pp.CFD_DEM["quadrature_order"].GetInt())

    if pp.CFD_DEM["drag_force_type"].GetInt() in {13} or pp.CFD_DEM["lift_force_type"].GetInt() == 1:
        dem_model_part.ProcessInfo.SetValue(POWER_LAW_K, pp.CFD_DEM["power_law_k"].GetDouble())
        dem_model_part.ProcessInfo.SetValue(POWER_LAW_N, pp.CFD_DEM["power_law_n"].GetDouble())


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
    pp.fluid_vars += [TORQUE]
    pp.fluid_vars += pp.fluid_printing_vars
    pp.fluid_vars += pp.coupling_fluid_vars

    if pp.CFD_DEM["pressure_grad_recovery_type"].GetInt() > 0:
        pp.fluid_vars += [RECOVERED_PRESSURE_GRADIENT]

    if pp.CFD_DEM["gradient_calculation_type"].GetInt() > 1 or pp.CFD_DEM["pressure_grad_recovery_type"].GetInt() > 1 or pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 7:
        pp.fluid_vars += [NODAL_WEIGHTS]

    if pp.CFD_DEM["material_acceleration_calculation_type"].GetInt():
        pp.fluid_vars += [MATERIAL_ACCELERATION]
        pp.fluid_vars += [VELOCITY_COMPONENT_GRADIENT]

        if pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 5 or pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 6:
            if pp.CFD_DEM["store_full_gradient_option"].GetBool():
                pp.fluid_vars += [VELOCITY_X_GRADIENT]
                pp.fluid_vars += [VELOCITY_Y_GRADIENT]
                pp.fluid_vars += [VELOCITY_Z_GRADIENT]

    if pp.CFD_DEM["vorticity_calculation_type"].GetInt() == 1 or pp.CFD_DEM["lift_force_type"].GetInt() == 1:
        pp.fluid_vars += [VORTICITY]

    if pp.CFD_DEM["laplacian_calculation_type"].GetInt():
        pp.fluid_vars += [VELOCITY_LAPLACIAN]

        if pp.CFD_DEM["include_faxen_terms_option"].GetBool():
            pp.fluid_vars += [VELOCITY_LAPLACIAN_RATE]

    if pp.CFD_DEM["calculate_diffusivity_option"].GetBool():
        pp.fluid_vars += [CONDUCTIVITY]

    # dem variables
    pp.dem_vars = []
    pp.dem_vars += pp.dem_printing_vars
    pp.dem_vars += pp.coupling_dem_vars
    pp.dem_vars += [BUOYANCY]
    pp.dem_vars += [VELOCITY_OLD]

    if pp.CFD_DEM["frame_of_reference_type"].GetInt() and pp.CFD_DEM["basset_force_type"].GetInt() > 0:
        pp.dem_vars += [DISPLACEMENT_OLD]
        pp.dem_vars += [VELOCITY_OLD_OLD]

    if pp.CFD_DEM["TranslationalIntegrationScheme"].GetString() in {'Hybrid_Bashforth', 'TerminalVelocityScheme'} or pp.CFD_DEM["basset_force_type"].GetInt() > 0:
        pp.dem_vars += [VELOCITY_OLD]
        pp.dem_vars += [ADDITIONAL_FORCE_OLD]
        pp.dem_vars += [SLIP_VELOCITY]

    if pp.CFD_DEM["drag_force_type"].GetInt() > 0 and  pp.CFD_DEM["add_each_hydro_force_option"].GetBool():
        pp.dem_vars += [DRAG_FORCE]

    if pp.CFD_DEM["drag_force_type"].GetInt() > 1:
        pp.dem_vars += [PARTICLE_SPHERICITY]

    if pp.CFD_DEM["lift_force_type"].GetInt() > 0 and  pp.CFD_DEM["add_each_hydro_force_option"].GetBool():
        pp.dem_vars += [LIFT_FORCE]

    if pp.CFD_DEM["virtual_mass_force_type"].GetInt() > 0 and  pp.CFD_DEM["add_each_hydro_force_option"].GetBool():
        pp.dem_vars += [VIRTUAL_MASS_FORCE]

    if pp.CFD_DEM["basset_force_type"].GetInt() > 0 and  pp.CFD_DEM["add_each_hydro_force_option"].GetBool():
        pp.dem_vars += [BASSET_FORCE]

    # clusters variables
    pp.clusters_vars = []

    # rigid faces variables
    pp.rigid_faces_vars = [VELOCITY, ANGULAR_VELOCITY, DISPLACEMENT, DELTA_DISPLACEMENT, DELTA_ROTATION, CONTACT_FORCES, DEM_PRESSURE, ELASTIC_FORCES, PRESSURE, TANGENTIAL_ELASTIC_FORCES, SHEAR_STRESS, NODAL_AREA]

    if pp.CFD_DEM["embedded_option"].GetBool():
        pp.rigid_faces_vars += [FORCE]
        pp.rigid_faces_vars += [POSITIVE_FACE_PRESSURE]
        pp.rigid_faces_vars += [NEGATIVE_FACE_PRESSURE]

    pp.fluid_vars += pp.rigid_faces_vars

    # inlet variables
    pp.inlet_vars = pp.dem_vars

def ConstructListsOfResultsToPrint(pp):
    pp.dem_nodal_results = []
    pp.clusters_nodal_results = []
    pp.rigid_faces_nodal_results = []

    if pp.CFD_DEM["print_SLIP_VELOCITY_option"].GetBool():
        pp.dem_nodal_results += ["SLIP_VELOCITY"]

    if pp.CFD_DEM["PostRadius"].GetBool():
        pp.dem_nodal_results += ["RADIUS"]

    if pp.CFD_DEM["PostAngularVelocity"].GetBool():
        pp.dem_nodal_results += ["ANGULAR_VELOCITY"]

    if pp.CFD_DEM["PostElasticForces"].GetBool():
        pp.dem_nodal_results += ["ELASTIC_FORCES"]

    if pp.CFD_DEM["PostContactForces"].GetBool():
        pp.dem_nodal_results += ["CONTACT_FORCES"]

    if pp.CFD_DEM["PostTotalForces"].GetBool():
        pp.dem_nodal_results += ["TOTAL_FORCES"]

    if pp.CFD_DEM["ElementType"].GetString() == "SwimmingNanoParticle":
        pp.dem_nodal_results += ["EXTERNAL_APPLIED_FORCE"]
        if pp.CFD_DEM.PostCationConcentration:
            pp.dem_nodal_results += ["CATION_CONCENTRATION"]

    if pp.CFD_DEM["coupling_level_type"].GetInt() > 0:

        if pp.CFD_DEM["print_REYNOLDS_NUMBER_option"].GetBool():
            pp.dem_nodal_results += ["REYNOLDS_NUMBER"]

        if pp.CFD_DEM["print_PRESSURE_GRAD_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["PRESSURE_GRAD_PROJECTED"]

        if pp.CFD_DEM["print_HYDRODYNAMIC_FORCE_option"].GetBool():
            pp.dem_nodal_results += ["HYDRODYNAMIC_FORCE"]

        if pp.CFD_DEM["print_HYDRODYNAMIC_MOMENT_option"].GetBool():
            pp.dem_nodal_results += ["HYDRODYNAMIC_MOMENT"]

        if pp.CFD_DEM["print_FLUID_VEL_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_VEL_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_VEL_PROJECTED_RATE_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_VEL_PROJECTED_RATE"]

        if pp.CFD_DEM["print_FLUID_VEL_LAPL_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_VEL_LAPL_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_VEL_LAPL_RATE_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_VEL_LAPL_RATE_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_ACCEL_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_ACCEL_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_FRACTION_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_FRACTION_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_FRACTION_GRADIENT_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_FRACTION_GRADIENT_PROJECTED"]

        if pp.CFD_DEM["print_FLUID_VISCOSITY_PROJECTED_option"].GetBool():
            pp.dem_nodal_results += ["FLUID_VISCOSITY_PROJECTED"]

        if pp.CFD_DEM["print_BUOYANCY_option"].GetBool():
            pp.dem_nodal_results += ["BUOYANCY"]

    if pp.CFD_DEM["add_each_hydro_force_option"].GetBool():

        if pp.CFD_DEM["print_DRAG_FORCE_option"].GetBool():
            pp.dem_nodal_results += ["DRAG_FORCE"]

        if pp.CFD_DEM["print_VIRTUAL_MASS_FORCE_option"].GetBool():
            pp.dem_nodal_results += ["VIRTUAL_MASS_FORCE"]

        if pp.CFD_DEM["print_BASSET_FORCE_option"].GetBool():
            pp.dem_nodal_results += ["BASSET_FORCE"]

        if pp.CFD_DEM["print_LIFT_FORCE_option"].GetBool():
            pp.dem_nodal_results += ["LIFT_FORCE"]

    if pp.CFD_DEM["embedded_option"].GetBool():
        pp.rigid_faces_nodal_results += ["POSITIVE_FACE_PRESSURE"]
        pp.rigid_faces_nodal_results += ["NEGATIVE_FACE_PRESSURE"]

    if pp.CFD_DEM["PostNonDimensionalVolumeWear"].GetBool():
        pp.rigid_faces_nodal_results += ["IMPACT_WEAR"]
        pp.rigid_faces_nodal_results += ["NON_DIMENSIONAL_VOLUME_WEAR"]

    # changes on the fluid variables to print for the sake of consistency
    ChangeListOfFluidNodalResultsToPrint(pp)

    pp.mixed_nodal_results = ["VELOCITY", "DISPLACEMENT"]

    pp.variables_to_print_in_file = ["DRAG_FORCE", "LIFT_FORCE", "BUOYANCY", "VELOCITY"]

    pp.dem_printing_vars = []
    pp.clusters_printing_vars = []
    pp.fluid_printing_vars = []
    pp.rigid_faces_printing_vars = []
    pp.time_filtered_vars = []

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

    if not pp.CFD_DEM["print_PRESSURE_option"].GetBool():
        if "PRESSURE" in pp.nodal_results:
            pp.nodal_results.remove("PRESSURE")

    EliminateRepeatedValuesFromList(pp.nodal_results)
    EliminateRepeatedValuesFromList(pp.dem_nodal_results)
    EliminateRepeatedValuesFromList(pp.mixed_nodal_results)

def ConstructListsOfVariablesForCoupling(pp):

    # fluid coupling variables
    pp.coupling_fluid_vars = []
    pp.coupling_fluid_vars += [MATERIAL_ACCELERATION]

    pp.coupling_fluid_vars += [KratosGlobals.GetVariable( pp.CFD_DEM["body_force_per_unit_mass_variable_name"].GetString() )]

    if pp.CFD_DEM["fluid_model_type"].GetInt() == 0:
        pp.coupling_fluid_vars += [AVERAGED_FLUID_VELOCITY]

    if pp.CFD_DEM["fluid_model_type"].GetInt() == 0 or pp.CFD_DEM["coupling_level_type"].GetInt() > 1 or pp.CFD_DEM["drag_force_type"].GetInt() == 4:
        pp.coupling_fluid_vars += [FLUID_FRACTION]
        pp.coupling_fluid_vars += [FLUID_FRACTION_OLD]

        if pp.CFD_DEM["print_DISPERSE_FRACTION_option"].GetBool():
            pp.coupling_fluid_vars += [DISPERSE_FRACTION]

        if pp.CFD_DEM["filter_velocity_option"].GetBool():
            pp.coupling_fluid_vars += [PARTICLE_VEL_FILTERED]
            pp.coupling_fluid_vars += [TIME_AVERAGED_ARRAY_3]
            pp.coupling_fluid_vars += [PHASE_FRACTION]

    if pp.CFD_DEM["fluid_model_type"].GetInt() >= 1:
        pp.coupling_fluid_vars += [FLUID_FRACTION_GRADIENT]
        pp.coupling_fluid_vars += [FLUID_FRACTION_RATE]

    if pp.CFD_DEM["coupling_level_type"].GetInt() >= 1:
        pp.coupling_fluid_vars += [HYDRODYNAMIC_REACTION]

    if pp.CFD_DEM["coupling_level_type"].GetInt() >= 1 and pp.CFD_DEM["time_averaging_type"].GetInt() > 0:
        pp.coupling_fluid_vars += [MEAN_HYDRODYNAMIC_REACTION]

    if pp.CFD_DEM["drag_force_type"].GetInt() in {2} or pp.CFD_DEM["lift_force_type"].GetInt() == 1:
        pp.coupling_fluid_vars += [POWER_LAW_N]
        pp.coupling_fluid_vars += [POWER_LAW_K]
        pp.coupling_fluid_vars += [YIELD_STRESS]

        if pp.CFD_DEM["drag_force_type"].GetInt() == 2:
            pp.coupling_fluid_vars += [GEL_STRENGTH]

    if pp.viscosity_modification_type:
        pp.coupling_fluid_vars += [VISCOSITY]

    if pp.CFD_DEM["embedded_option"].GetBool():
        pp.coupling_fluid_vars += [DISTANCE]

    # dem coupling variables
    pp.coupling_dem_vars = []

    if pp.CFD_DEM["coupling_level_type"].GetInt() > 0:
        pp.coupling_dem_vars += [FLUID_VEL_PROJECTED]
        pp.coupling_dem_vars += [FLUID_DENSITY_PROJECTED]
        pp.coupling_dem_vars += [FLUID_VISCOSITY_PROJECTED]
        pp.coupling_dem_vars += [HYDRODYNAMIC_FORCE]
        pp.coupling_dem_vars += [HYDRODYNAMIC_MOMENT]
        pp.coupling_dem_vars += [MATERIAL_FLUID_ACCEL_PROJECTED]
        pp.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]
        pp.coupling_dem_vars += [FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED]
        pp.coupling_dem_vars += [ADDITIONAL_FORCE] # Here for safety for the moment

        if pp.CFD_DEM["buoyancy_force_type"].GetInt() != 2 and pp.CFD_DEM["drag_force_type"].GetInt() != 2:
            pp.coupling_dem_vars += [PRESSURE_GRAD_PROJECTED]

        if pp.CFD_DEM["include_faxen_terms_option"].GetBool():
            pp.coupling_dem_vars += [FLUID_VEL_LAPL_PROJECTED]
            pp.coupling_dem_vars += [FLUID_VEL_LAPL_RATE_PROJECTED]

        if pp.CFD_DEM["basset_force_type"].GetInt() > 0:
            pp.coupling_dem_vars += [FLUID_VEL_PROJECTED_RATE]

    if pp.CFD_DEM["coupling_level_type"].GetInt() >= 1 or pp.CFD_DEM["fluid_model_type"].GetInt() == 0:
        pp.coupling_dem_vars += [FLUID_FRACTION_PROJECTED]

    if pp.CFD_DEM["coupling_level_type"].GetInt() >= 1 and pp.CFD_DEM["print_FLUID_FRACTION_GRADIENT_PROJECTED_option"].GetBool():
        pp.coupling_dem_vars += [FLUID_FRACTION_GRADIENT_PROJECTED]

    if pp.CFD_DEM["drag_force_type"].GetInt() in {2} or pp.CFD_DEM["lift_force_type"].GetInt() == 1:
        pp.coupling_dem_vars += [POWER_LAW_N]
        pp.coupling_dem_vars += [POWER_LAW_K]
        pp.coupling_dem_vars += [YIELD_STRESS]

    if pp.CFD_DEM["lift_force_type"].GetInt() == 1:
        pp.coupling_dem_vars += [FLUID_VORTICITY_PROJECTED]
        pp.coupling_dem_vars += [SHEAR_RATE_PROJECTED]

    if pp.CFD_DEM["virtual_mass_force_type"].GetInt() or pp.CFD_DEM["basset_force_type"].GetInt() >= 1:
        pp.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]

    if pp.CFD_DEM["embedded_option"].GetBool():
        pp.coupling_dem_vars += [DISTANCE]

    if pp.CFD_DEM["print_REYNOLDS_NUMBER_option"].GetBool():
        pp.coupling_dem_vars += [REYNOLDS_NUMBER]

    if pp.CFD_DEM["apply_time_filter_to_fluid_fraction_option"].GetBool():
        pp.time_filtered_vars += [FLUID_FRACTION_FILTERED]

    if pp.CFD_DEM["filter_velocity_option"].GetBool():
        pp.time_filtered_vars += [PARTICLE_VEL_FILTERED]


def ChangeListOfFluidNodalResultsToPrint(pp):

    if pp.CFD_DEM["store_full_gradient_option"].GetBool() and pp.CFD_DEM["print_VELOCITY_GRADIENT_option"].GetBool():
        pp.nodal_results += ["VELOCITY_X_GRADIENT"]
        pp.nodal_results += ["VELOCITY_Y_GRADIENT"]
        pp.nodal_results += ["VELOCITY_Z_GRADIENT"]

    if pp.CFD_DEM["print_PRESSURE_GRADIENT_option"].GetBool():
        pp.nodal_results += ["PRESSURE_GRADIENT"]

    if pp.CFD_DEM["print_FLUID_FRACTION_option"].GetBool():
        pp.nodal_results += ["FLUID_FRACTION"]

    if pp.CFD_DEM["print_PARTICLE_VEL_option"].GetBool():
        pp.nodal_results += ["PARTICLE_VEL_FILTERED"]

    if pp.CFD_DEM["print_FLUID_FRACTION_GRADIENT_option"].GetBool():
        pp.nodal_results += ["FLUID_FRACTION_GRADIENT"]

    if pp.CFD_DEM["print_DISPERSE_FRACTION_option"].GetBool():
        pp.nodal_results += ["DISPERSE_FRACTION"]
        pp.nodal_results += ["PHASE_FRACTION"]

    if pp.CFD_DEM["fluid_model_type"].GetInt() == 0 and pp.CFD_DEM["print_AVERAGED_FLUID_VELOCITY_option"].GetBool():
        pp.nodal_results += ["AVERAGED_FLUID_VELOCITY"]

    if pp.CFD_DEM["fluid_model_type"].GetInt() == 1 and pp.CFD_DEM["print_FLUID_FRACTION_GRADIENT_option"].GetBool():
        pp.nodal_results += ["FLUID_FRACTION_GRADIENT"]

    if pp.CFD_DEM["print_VISCOSITY_option"].GetBool():
        pp.nodal_results += ["VISCOSITY"]

    if pp.CFD_DEM["body_force_on_fluid_option"].GetBool() and pp.CFD_DEM["print_BODY_FORCE_option"].GetBool():
        pp.nodal_results += ["BODY_FORCE"]

    if pp.CFD_DEM["print_HYDRODYNAMIC_REACTION_option"].GetBool():
        pp.nodal_results += ["HYDRODYNAMIC_REACTION"]

    if pp.CFD_DEM["print_MEAN_HYDRODYNAMIC_REACTION_option"].GetBool():
        pp.nodal_results += ["MEAN_HYDRODYNAMIC_REACTION"]

    if pp.CFD_DEM["embedded_option"].GetBool() and pp.CFD_DEM["print_distance_option"].GetBool():
        pp.nodal_results += ["DISTANCE"]

    if pp.CFD_DEM["print_MATERIAL_ACCELERATION_option"].GetBool():
        pp.nodal_results += ["MATERIAL_ACCELERATION"]

    if pp.CFD_DEM["print_VORTICITY_option"].GetBool():
        pp.nodal_results += ["VORTICITY"]

    if pp.CFD_DEM["print_VELOCITY_LAPLACIAN_option"].GetBool():
        pp.nodal_results += ["VELOCITY_LAPLACIAN"]

    if pp.CFD_DEM["print_VELOCITY_LAPLACIAN_RATE_option"].GetBool():
        pp.nodal_results += ["VELOCITY_LAPLACIAN_RATE"]

    if pp.CFD_DEM["print_CONDUCTIVITY_option"].GetBool():
        pp.nodal_results += ["CONDUCTIVITY"]

    if pp.CFD_DEM["print_VECTORIAL_ERROR_option"].GetBool():
        pp.nodal_results += ["VECTORIAL_ERROR"]
        pp.nodal_results += ["VECTORIAL_ERROR_1"]

def ChangeInputDataForConsistency(pp):
    if pp.CFD_DEM["coupling_level_type"].GetInt() == 0:
        pp.CFD_DEM["project_at_every_substep_option"].SetBool(False)
        pp.CFD_DEM.project_at_every_substep_option = False #TODO: what should we do in these cases?? Discuss

    if pp.CFD_DEM["VelocityTrapOption"]:
        pp.velocity_trap_option = 1
    else:
        pp.velocity_trap_option = 0

    if pp.CFD_DEM["flow_in_porous_medium_option"].GetBool():
        pp.coupling_weighing_type = - 1 # the fluid fraction is not projected from DEM (there may not be a DEM part) but is externally imposed

    pp.CFD_DEM["time_steps_per_stationarity_step"].SetInt(max(1, int(pp.CFD_DEM["time_steps_per_stationarity_step"].GetInt())))

    if pp.CFD_DEM["coupling_level_type"].GetInt() > 1:
        pp.CFD_DEM["stationary_problem_option"].SetBool(False)

def EliminateRepeatedValuesFromList(redundant_list):
    clean_list = []

    for var in redundant_list:

        if var in clean_list:
            redundant_list.remove(var)

        clean_list += [var]
