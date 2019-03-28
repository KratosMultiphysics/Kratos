from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import ast
import parameters_tools as PT
def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)

class VariablesManager:
    @staticmethod
    def EliminateRepeatedValuesFromList(redundant_list):
        clean_list = []

        for var in redundant_list:

            if var in clean_list:
                redundant_list.remove(var)

            clean_list += [var]

    @staticmethod
    def AddNodalVariables(model_part, variable_list):

        for var in variable_list:
            model_part.AddNodalSolutionStepVariable(var)

    def __init__(self, parameters):
        self.project_parameters = parameters

    # constructing lists of variables to add
    # * Performing modifications to the input parameters for consistency (provisional until interface does it)
    # * Choosing the variables to be printed
    # * Choosing the variables to be passed as a parameter to the constructor of a ProjectionModule
    #       instance to be filled with the other phase's info through the coupling process
    # * Listing nodal variables to be added to the model parts (memory will be allocated for them).
    #       Note that additional variables may be added as well by the fluid and/or DEM strategies.
    @staticmethod
    def AddFrameOfReferenceRelatedVariables(parameters, model_part):
        frame_of_reference_type = parameters["frame_of_reference_type"].GetInt()
        model_part.ProcessInfo.SetValue(FRAME_OF_REFERENCE_TYPE, frame_of_reference_type)

        if frame_of_reference_type == 1: # Rotating frame
            angular_velocity_of_frame = Vector(3)
            angular_velocity_of_frame[:] = [parameters["angular_velocity_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]

            model_part.ProcessInfo.SetValue(ANGULAR_VELOCITY_MOVING_FRAME, angular_velocity_of_frame)

            if frame_of_reference_type >= 2: # Gemeral frame
                angular_velocity_of_frame_old = Vector(3)
                angular_velocity_of_frame_old[:] = [parameters["angular_velocity_of_frame_old" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                acceleration_of_frame_origin = Vector(3)
                acceleration_of_frame_origin[:] = [parameters["acceleration_of_frame_origin" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                angular_acceleration_of_frame = Vector(3)
                angular_acceleration_of_frame[:] = [parameters["angular_acceleration_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                model_part.ProcessInfo.SetValue(ANGULAR_VELOCITY_MOVING_FRAME_OLD, angular_velocity_of_frame_old)
                model_part.ProcessInfo.SetValue(ACCELERATION_MOVING_FRAME_ORIGIN, acceleration_of_frame_origin)
                model_part.ProcessInfo.SetValue(ANGULAR_ACCELERATION_MOVING_FRAME, angular_acceleration_of_frame)

    def AddExtraProcessInfoVariablesToFluidModelPart(self, parameters, fluid_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, fluid_model_part)

        fluid_model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1)
        gravity = Vector(3)
        if parameters["body_force_on_fluid_option"].GetBool():
            gravity[0] = parameters["GravityX"].GetDouble()
            gravity[1] = parameters["GravityY"].GetDouble()
            gravity[2] = parameters["GravityZ"].GetDouble()
        fluid_model_part.ProcessInfo.SetValue(GRAVITY, gravity)

        if parameters["laplacian_calculation_type"].GetInt() == 3: # recovery through solving a system
            fluid_model_part.ProcessInfo.SetValue(COMPUTE_LUMPED_MASS_MATRIX, 1)
        elif (parameters["material_acceleration_calculation_type"].GetInt() == 4
              or parameters["material_acceleration_calculation_type"].GetInt() == 5
              or parameters["material_acceleration_calculation_type"].GetInt() == 6): # recovery by solving a system
            fluid_model_part.ProcessInfo.SetValue(COMPUTE_LUMPED_MASS_MATRIX, 0)

        if parameters["material_acceleration_calculation_type"].GetInt() == 5 or parameters["material_acceleration_calculation_type"].GetInt() == 6:
            fluid_model_part.ProcessInfo.SetValue(CURRENT_COMPONENT, 0)

        if parameters["non_newtonian_option"].GetBool():
            fluid_model_part.ProcessInfo.SetValue(YIELD_STRESS, parameters["yield_stress"].GetDouble())
            fluid_model_part.ProcessInfo.SetValue(REGULARIZATION_COEFFICIENT, parameters["regularization_coefficient"].GetDouble())
            fluid_model_part.ProcessInfo.SetValue(POWER_LAW_K, parameters["power_law_k"].GetDouble())
            fluid_model_part.ProcessInfo.SetValue(POWER_LAW_N, parameters["power_law_n"].GetDouble())

    def AddExtraProcessInfoVariablesToDispersePhaseModelPart(self, parameters, dem_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, dem_model_part)
        dem_model_part.ProcessInfo.SetValue(COUPLING_TYPE, parameters["coupling_level_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(BASSET_FORCE_TYPE, parameters["basset_force_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, parameters["fluid_model_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, parameters["drag_modifier_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, parameters["initial_drag_force"].GetDouble())
        dem_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, parameters["drag_law_slope"].GetDouble())
        dem_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, parameters["power_law_tol"].GetDouble())
        dem_model_part.ProcessInfo.SetValue(DRAG_POROSITY_CORRECTION_TYPE, parameters["drag_porosity_correction_type"].GetInt())

        for prop in parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]

                if (prop["hydrodynamic_law_parameters"]["name"].GetString() != "default"
                    and history_force_parameters["name"].GetString() != "default"):
                    dem_model_part.ProcessInfo.SetValue(NUMBER_OF_INIT_BASSET_STEPS, history_force_parameters["n_init_basset_steps"].GetInt())
                    dem_model_part.ProcessInfo.SetValue(TIME_STEPS_PER_QUADRATURE_STEP, history_force_parameters["time_steps_per_quadrature_step"].GetInt())
                    dem_model_part.ProcessInfo.SetValue(LAST_TIME_APPENDING, 0.0)
                    dem_model_part.ProcessInfo.SetValue(QUADRATURE_ORDER, history_force_parameters["quadrature_order"].GetInt())
                    break

        if parameters["non_newtonian_option"].GetBool():
            dem_model_part.ProcessInfo.SetValue(POWER_LAW_K, parameters["power_law_k"].GetDouble())
            dem_model_part.ProcessInfo.SetValue(POWER_LAW_N, parameters["power_law_n"].GetDouble())

    def ConstructListsOfVariables(self, parameters):
        # PRINTING VARIABLES
        # constructing lists of variables to be printed
        self.ConstructListsOfResultsToPrint(parameters)

        # COUPLING VARIABLES
        # listing the variables involved in the fluid-particles coupling

        if parameters["coupling_level_type"].GetInt():
            self.ConstructListsOfVariablesForCoupling(parameters)

        # VARIABLES TO ADD
        # listing nodal variables to be added to the model parts (memory will be allocated for them)

        # fluid variables
        self.fluid_vars = []
        self.fluid_vars += [TORQUE]
        self.fluid_vars += self.fluid_printing_vars
        self.fluid_vars += self.coupling_fluid_vars

        if parameters["pressure_grad_recovery_type"].GetInt() > 0:
            self.fluid_vars += [RECOVERED_PRESSURE_GRADIENT]

        if (parameters["gradient_calculation_type"].GetInt() > 1
            or parameters["pressure_grad_recovery_type"].GetInt() > 1
            or parameters["material_acceleration_calculation_type"].GetInt() == 7
            or parameters["laplacian_calculation_type"].GetInt() > 1):
            self.fluid_vars += [NODAL_WEIGHTS]

        if parameters["material_acceleration_calculation_type"].GetInt():
            self.fluid_vars += [MATERIAL_ACCELERATION]
            self.fluid_vars += [VELOCITY_COMPONENT_GRADIENT]

            if (parameters["material_acceleration_calculation_type"].GetInt() == 5
                or parameters["material_acceleration_calculation_type"].GetInt() == 6):
                if parameters["store_full_gradient_option"].GetBool():
                    self.fluid_vars += [VELOCITY_X_GRADIENT]
                    self.fluid_vars += [VELOCITY_Y_GRADIENT]
                    self.fluid_vars += [VELOCITY_Z_GRADIENT]

        if (parameters["vorticity_calculation_type"].GetInt() > 0
            or PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')):
            self.fluid_vars += [VORTICITY]

        if parameters["laplacian_calculation_type"].GetInt():
            self.fluid_vars += [VELOCITY_LAPLACIAN]

        if PT.RecursiveFindTrueBoolInParameters(parameters["properties"], 'do_apply_faxen_corrections'):
            self.fluid_vars += [VELOCITY_LAPLACIAN_RATE]

        if parameters["calculate_diffusivity_option"].GetBool():
            self.fluid_vars += [CONDUCTIVITY]

        # dem variables
        self.dem_vars = []
        self.dem_vars += self.dem_printing_vars
        self.dem_vars += self.coupling_dem_vars
        self.dem_vars += [BUOYANCY]
        self.dem_vars += [VELOCITY_OLD]

        if parameters["frame_of_reference_type"].GetInt() and parameters["basset_force_type"].GetInt() > 0:
            self.dem_vars += [DISPLACEMENT_OLD]
            self.dem_vars += [VELOCITY_OLD_OLD]

        if (parameters["TranslationalIntegrationScheme"].GetString()
            in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}
            or parameters["basset_force_type"].GetInt() > 0):
            self.dem_vars += [VELOCITY_OLD]
            self.dem_vars += [ADDITIONAL_FORCE_OLD]
            self.dem_vars += [AUX_VEL]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [DRAG_FORCE]

        self.dem_vars += [PARTICLE_SPHERICITY] # TODO: add only when needed

        if (PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')
            and parameters["add_each_hydro_force_option"].GetBool()):
            self.dem_vars += [LIFT_FORCE]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [VIRTUAL_MASS_FORCE]

        will_need_basset_force_variable = False
        for prop in parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):
                if prop["hydrodynamic_law_parameters"]["history_force_parameters"]["name"].GetString() != 'default':
                    will_need_basset_force_variable = True
                    break

        if will_need_basset_force_variable:
            self.dem_vars += [BASSET_FORCE]

        # clusters variables
        self.clusters_vars = []

        # rigid faces variables
        self.rigid_faces_vars = [VELOCITY,
                                 ANGULAR_VELOCITY,
                                 DISPLACEMENT,
                                 DELTA_DISPLACEMENT,
                                 DELTA_ROTATION,
                                 CONTACT_FORCES,
                                 DEM_PRESSURE,
                                 ELASTIC_FORCES,
                                 PRESSURE,
                                 TANGENTIAL_ELASTIC_FORCES,
                                 SHEAR_STRESS,
                                 NODAL_AREA,
                                 VELOCITY_OLD]

        if parameters["embedded_option"].GetBool():
            self.rigid_faces_vars += [FORCE]
            self.rigid_faces_vars += [POSITIVE_FACE_PRESSURE]
            self.rigid_faces_vars += [NEGATIVE_FACE_PRESSURE]

        self.fluid_vars += self.rigid_faces_vars

        # inlet variables
        self.inlet_vars = self.dem_vars

    def ConstructListsOfResultsToPrint(self, parameters):
        dem_list = self.project_parameters["dem_nodal_results"]
        self.dem_nodal_results = [key for key in dem_list.keys() if dem_list[key].GetBool()]
        self.clusters_nodal_results = []
        self.rigid_faces_nodal_results = []

        if parameters['dem_parameters']["PostRadius"].GetBool():
            self.dem_nodal_results += ["RADIUS"]

        if parameters['dem_parameters']["PostAngularVelocity"].GetBool():
            self.dem_nodal_results += ["ANGULAR_VELOCITY"]

        if parameters['dem_parameters']["PostElasticForces"].GetBool():
            self.dem_nodal_results += ["ELASTIC_FORCES"]

        if parameters['dem_parameters']["PostContactForces"].GetBool():
            self.dem_nodal_results += ["CONTACT_FORCES"]

        if parameters['dem_parameters']["PostTotalForces"].GetBool():
            self.dem_nodal_results += ["TOTAL_FORCES"]

        if parameters["ElementType"].GetString() == "SwimmingNanoParticle":
            self.dem_nodal_results += ["EXTERNAL_APPLIED_FORCE"]
            if parameters.PostCationConcentration:
                self.dem_nodal_results += ["CATION_CONCENTRATION"]

        if parameters["embedded_option"].GetBool():
            self.rigid_faces_nodal_results += ["POSITIVE_FACE_PRESSURE"]
            self.rigid_faces_nodal_results += ["NEGATIVE_FACE_PRESSURE"]

        if parameters['dem_parameters']["PostNonDimensionalVolumeWear"].GetBool():
            self.rigid_faces_nodal_results += ["IMPACT_WEAR"]
            self.rigid_faces_nodal_results += ["NON_DIMENSIONAL_VOLUME_WEAR"]

        # changes on the fluid variables to print for the sake of consistency
        self.ChangeListOfFluidNodalResultsToPrint(parameters)

        self.mixed_nodal_results = ["VELOCITY", "DISPLACEMENT"]

        self.variables_to_print_in_file = ["DRAG_FORCE", "LIFT_FORCE", "BUOYANCY", "VELOCITY"]

        self.dem_printing_vars = []

        self.clusters_printing_vars = []
        self.fluid_printing_vars = []
        self.rigid_faces_printing_vars = []
        self.time_filtered_vars = []

        for variable in self.nodal_results:
            self.fluid_printing_vars += [eval(variable)]

        for variable in self.dem_nodal_results:
            self.dem_printing_vars += [eval(variable)]

        for variable in self.clusters_nodal_results:
            self.clusters_printing_vars += [eval(variable)]

        for variable in self.rigid_faces_nodal_results:
            self.rigid_faces_printing_vars += [eval(variable)]

        for variable in self.mixed_nodal_results:
            self.dem_printing_vars += [eval(variable)]
            self.fluid_printing_vars += [eval(variable)]

        for var in self.mixed_nodal_results:

            if var in self.nodal_results:
                self.nodal_results.remove(var)

        VariablesManager.EliminateRepeatedValuesFromList(self.nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.dem_nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.mixed_nodal_results)

    def ConstructListsOfVariablesForCoupling(self, parameters):

        # fluid coupling variables
        self.coupling_fluid_vars = []
        self.coupling_fluid_vars += [MATERIAL_ACCELERATION]

        self.coupling_fluid_vars += [KratosGlobals.GetVariable(parameters["body_force_per_unit_mass_variable_name"].GetString() )]

        if parameters["fluid_model_type"].GetInt() == 0:
            self.coupling_fluid_vars += [AVERAGED_FLUID_VELOCITY]

        if (parameters["fluid_model_type"].GetInt() == 0
            or parameters["coupling_level_type"].GetInt() > 1):
            self.coupling_fluid_vars += [FLUID_FRACTION]
            self.coupling_fluid_vars += [FLUID_FRACTION_OLD]

            if 'DISPERSE_FRACTION' in self.nodal_results:
                self.coupling_fluid_vars += [DISPERSE_FRACTION]

            if parameters["filter_velocity_option"].GetBool():
                self.coupling_fluid_vars += [PARTICLE_VEL_FILTERED]
                self.coupling_fluid_vars += [TIME_AVERAGED_ARRAY_3]
                self.coupling_fluid_vars += [PHASE_FRACTION]

        if parameters["fluid_model_type"].GetInt() >= 1:
            self.coupling_fluid_vars += [FLUID_FRACTION_GRADIENT]
            self.coupling_fluid_vars += [FLUID_FRACTION_RATE]

        if parameters["coupling_level_type"].GetInt() >= 1:
            self.coupling_fluid_vars += [HYDRODYNAMIC_REACTION]

        if parameters["coupling_level_type"].GetInt() >= 1 and parameters["time_averaging_type"].GetInt() > 0:
            self.coupling_fluid_vars += [MEAN_HYDRODYNAMIC_REACTION]

        if parameters["non_newtonian_option"].GetBool():
            self.coupling_fluid_vars += [POWER_LAW_N]
            self.coupling_fluid_vars += [POWER_LAW_K]
            self.coupling_fluid_vars += [YIELD_STRESS]
            # self.coupling_fluid_vars += [GEL_STRENGTH] # TODO: make specific option for this

        if parameters["viscosity_modification_type"].GetInt():
            self.coupling_fluid_vars += [VISCOSITY]

        if parameters["embedded_option"].GetBool():
            self.coupling_fluid_vars += [DISTANCE]

        # dem coupling variables
        self.coupling_dem_vars = []

        if parameters["coupling_level_type"].GetInt() > 0:
            self.coupling_dem_vars += [FLUID_VEL_PROJECTED]
            self.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]
            self.coupling_dem_vars += [FLUID_DENSITY_PROJECTED]
            self.coupling_dem_vars += [FLUID_VISCOSITY_PROJECTED]
            self.coupling_dem_vars += [HYDRODYNAMIC_FORCE]
            self.coupling_dem_vars += [HYDRODYNAMIC_MOMENT]
            self.coupling_dem_vars += [MATERIAL_FLUID_ACCEL_PROJECTED]
            self.coupling_dem_vars += [FLUID_ACCEL_PROJECTED]
            self.coupling_dem_vars += [FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED]
            self.coupling_dem_vars += [ADDITIONAL_FORCE] # Here for safety for the moment

            if PT.RecursiveFindTrueBoolInParameters(parameters["properties"], 'do_apply_faxen_corrections'):
                self.coupling_dem_vars += [FLUID_VEL_LAPL_PROJECTED]
                self.coupling_dem_vars += [FLUID_VEL_LAPL_RATE_PROJECTED]

            if parameters["basset_force_type"].GetInt() > 0:
                self.coupling_dem_vars += [FLUID_VEL_PROJECTED_RATE]

        if parameters["coupling_level_type"].GetInt() >= 1 or parameters["fluid_model_type"].GetInt() == 0:
            self.coupling_dem_vars += [FLUID_FRACTION_PROJECTED]

        if (parameters["coupling_level_type"].GetInt() >= 1
            and 'FLUID_FRACTION_GRADIENT_PROJECTED' in self.dem_printing_vars):
            self.coupling_dem_vars += [FLUID_FRACTION_GRADIENT_PROJECTED]

        if parameters["non_newtonian_option"].GetBool():
            self.coupling_dem_vars += [POWER_LAW_N]
            self.coupling_dem_vars += [POWER_LAW_K]
            self.coupling_dem_vars += [YIELD_STRESS]

        if PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters'):
            self.coupling_dem_vars += [FLUID_VORTICITY_PROJECTED]
            self.coupling_dem_vars += [SHEAR_RATE_PROJECTED]

        if parameters["embedded_option"].GetBool():
            self.coupling_dem_vars += [DISTANCE]

        if 'REYNOLDS_NUMBER' in self.dem_nodal_results:
            self.coupling_dem_vars += [REYNOLDS_NUMBER]

        if parameters["apply_time_filter_to_fluid_fraction_option"].GetBool():
            self.time_filtered_vars += [FLUID_FRACTION_FILTERED]

        if parameters["filter_velocity_option"].GetBool():
            self.time_filtered_vars += [PARTICLE_VEL_FILTERED]


    def ChangeListOfFluidNodalResultsToPrint(self, parameters):

        if parameters["store_full_gradient_option"].GetBool() and 'VELOCITY_GRADIENT' in self.nodal_results:
            self.nodal_results += ["VELOCITY_X_GRADIENT"]
            self.nodal_results += ["VELOCITY_Y_GRADIENT"]
            self.nodal_results += ["VELOCITY_Z_GRADIENT"]

        if parameters["fluid_model_type"].GetInt() == 0 and 'AVERAGED_FLUID_VELOCITY' in self.nodal_results:
            self.nodal_results += ["AVERAGED_FLUID_VELOCITY"]

        if parameters["fluid_model_type"].GetInt() == 1 and 'FLUID_FRACTION_GRADIENT' in self.nodal_results:
            self.nodal_results += ["FLUID_FRACTION_GRADIENT"]