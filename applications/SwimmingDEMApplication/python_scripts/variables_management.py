import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, Logger
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.FluidDynamicsApplication as Fluid
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT

def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)

def GetGlobalVariableByName(variable_name):
    modules = [Kratos, DEM, SDEM]
    for mod in modules:
        try:
            return Kratos.KratosGlobals.GetVariable(variable_name)
        except Exception:
            pass
    names = [mod.__name__ for mod in modules]
    error_message = ('No variable with name \''
        + variable_name + '\' exists in either of the modules:\n')
    for name in names[:-1]:
        error_message += name + ', '
    error_message += 'or ' + names[-1] + '.'
    raise AttributeError(error_message)

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
        self.SetOptions(parameters)
    # constructing lists of variables to add
    # * Performing modifications to the input parameters for consistency (provisional until interface does it)
    # * Choosing the variables to be printed
    # * Choosing the variables to be passed as a parameter to the constructor of a ProjectionModule
    #       instance to be filled with the other phase's info through the coupling process
    # * Listing nodal variables to be added to the model parts (memory will be allocated for them).
    #       Note that additional variables may be added as well by the fluid and/or DEM strategies.
    @staticmethod
    def AddFrameOfReferenceRelatedVariables(parameters, model_part):
        frame_of_reference_type = parameters["frame_of_reference"]["frame_type"].GetInt()
        model_part.ProcessInfo.SetValue(Kratos.FRAME_OF_REFERENCE_TYPE, frame_of_reference_type)

        if frame_of_reference_type == 1: # Rotating frame
            angular_velocity_of_frame = Vector(3)
            angular_velocity_of_frame[:] = [parameters['frame_of_reference']["angular_velocity_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]

            model_part.ProcessInfo.SetValue(Kratos.ANGULAR_VELOCITY_MOVING_FRAME, angular_velocity_of_frame)

            if frame_of_reference_type >= 2: # Gemeral frame
                angular_velocity_of_frame_old = Vector(3)
                angular_velocity_of_frame_old[:] = [parameters['frame_of_reference']["angular_velocity_of_frame_old" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                acceleration_of_frame_origin = Vector(3)
                acceleration_of_frame_origin[:] = [parameters['frame_of_reference']["acceleration_of_frame_origin" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                angular_acceleration_of_frame = Vector(3)
                angular_acceleration_of_frame[:] = [parameters['frame_of_reference']["angular_acceleration_of_frame" + comp].GetDouble() for comp in ['_X', '_Y', '_Z']][:]
                model_part.ProcessInfo.SetValue(Kratos.ANGULAR_VELOCITY_MOVING_FRAME_OLD, angular_velocity_of_frame_old)
                model_part.ProcessInfo.SetValue(Kratos.ACCELERATION_MOVING_FRAME_ORIGIN, acceleration_of_frame_origin)
                model_part.ProcessInfo.SetValue(Kratos.ANGULAR_ACCELERATION_MOVING_FRAME, angular_acceleration_of_frame)


    def SetOptions(self, parameters):
        self.do_include_history_force = (PT.RecursiveFindParametersWithCondition(
                                         parameters["properties"], 'history_force_parameters',
                                         condition=lambda value: value['name'].GetString() != 'default'))

        if self.do_include_history_force: #TODO: extend to multiple properties
            for prop in parameters["properties"].values():
                self.history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                break

        self.do_backward_coupling = parameters["coupling"]["coupling_level_type"].GetInt() > 1
        self.backward_coupling_tools = parameters["coupling"]["backward_coupling"]

    def AddExtraProcessInfoVariablesToFluidModelPart(self, parameters, fluid_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, fluid_model_part)

        fluid_model_part.ProcessInfo.SetValue(Kratos.FRACTIONAL_STEP, 1)

        if parameters["custom_fluid"]["body_force_on_fluid_option"].GetBool():
            gravity = self.project_parameters["gravity_parameters"]["direction"].GetVector()
            gravity *= self.project_parameters["gravity_parameters"]["modulus"].GetDouble()
            fluid_model_part.ProcessInfo.SetValue(Kratos.GRAVITY, gravity)
        else:
            fluid_model_part.ProcessInfo.SetValue(Kratos.GRAVITY, Vector([0.0] * 3))

        if parameters["laplacian_calculation_type"].GetInt() == 3: # recovery through solving a system
            fluid_model_part.ProcessInfo.SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, 1)
        elif (parameters["material_acceleration_calculation_type"].GetInt() == 4
              or parameters["material_acceleration_calculation_type"].GetInt() == 5
              or parameters["material_acceleration_calculation_type"].GetInt() == 6): # recovery by solving a system
            fluid_model_part.ProcessInfo.SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, 0)

        if parameters["material_acceleration_calculation_type"].GetInt() == 5 or parameters["material_acceleration_calculation_type"].GetInt() == 6:
            fluid_model_part.ProcessInfo.SetValue(Kratos.CURRENT_COMPONENT, 0)

        #TODO: review this lines
        # if parameters["non_newtonian_fluid"]["non_newtonian_option"].GetBool():
        #     fluid_model_part.ProcessInfo.SetValue(Kratos.YIELD_STRESS, parameters["non_newtonian_fluid"]["yield_stress"].GetDouble())
        #     fluid_model_part.ProcessInfo.SetValue(Kratos.REGULARIZATION_COEFFICIENT, parameters["non_newtonian_fluid"]["regularization_coefficient"].GetDouble())
        #     fluid_model_part.ProcessInfo.SetValue(Kratos.POWER_LAW_K, parameters["non_newtonian_fluid"]["power_law_k"].GetDouble())
        #     fluid_model_part.ProcessInfo.SetValue(Kratos.POWER_LAW_N, parameters["non_newtonian_fluid"]["power_law_n"].GetDouble())

    def AddExtraProcessInfoVariablesToDispersePhaseModelPart(self, parameters, dem_model_part):

        VariablesManager.AddFrameOfReferenceRelatedVariables(parameters, dem_model_part)
        dem_model_part.ProcessInfo.SetValue(Kratos.COUPLING_TYPE, parameters["coupling"]["coupling_level_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(Kratos.FLUID_MODEL_TYPE, parameters["custom_fluid"]["fluid_model_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(Kratos.DRAG_MODIFIER_TYPE, parameters["drag_modifier_type"].GetInt())
        dem_model_part.ProcessInfo.SetValue(Kratos.POWER_LAW_TOLERANCE, parameters["non_newtonian_fluid"]["power_law_tol"].GetDouble())

        if self.do_include_history_force:
            dem_model_part.ProcessInfo.SetValue(Kratos.NUMBER_OF_INIT_BASSET_STEPS, self.history_force_parameters["n_init_basset_steps"].GetInt())
            dem_model_part.ProcessInfo.SetValue(Kratos.TIME_STEPS_PER_QUADRATURE_STEP, self.history_force_parameters["time_steps_per_quadrature_step"].GetInt())
            dem_model_part.ProcessInfo.SetValue(Kratos.LAST_TIME_APPENDING, 0.0)
            dem_model_part.ProcessInfo.SetValue(Kratos.QUADRATURE_ORDER, self.history_force_parameters["quadrature_order"].GetInt())

        if parameters["non_newtonian_fluid"]["non_newtonian_option"].GetBool():
            dem_model_part.ProcessInfo.SetValue(Kratos.POWER_LAW_K, parameters["non_newtonian_fluid"]["power_law_k"].GetDouble())
            dem_model_part.ProcessInfo.SetValue(Kratos.POWER_LAW_N, parameters["non_newtonian_fluid"]["power_law_n"].GetDouble())

    def ConstructListsOfVariables(self, parameters):
        # PRINTING VARIABLES
        # constructing lists of variables to be printed

        self.nodal_results, self.gauss_points_results = [], []
        self.fluid_parameters = parameters['fluid_parameters']
        if self.fluid_parameters.Has('sdem_output_processes'):
            gid_output_options = self.fluid_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            gauss_point_results = result_file_configuration["gauss_point_results"]
            nodal_variables = self.fluid_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"]
            self.nodal_results = [nodal_variables[i].GetString() for i in range(nodal_variables.size())]
            self.gauss_points_results = [gauss_point_results[i].GetString() for i in range(gauss_point_results.size())]

        self.ConstructListsOfResultsToPrint(parameters)

        # COUPLING VARIABLES
        # listing the variables involved in the fluid-particles coupling

        if parameters["coupling"]["coupling_level_type"].GetInt():
            self.ConstructListsOfVariablesForCoupling(parameters)

        # VARIABLES TO ADD
        # listing nodal variables to be added to the model parts (memory will be allocated for them)

        # fluid variables
        self.fluid_vars = []
        self.fluid_vars += [Kratos.TORQUE]
        self.fluid_vars += self.fluid_printing_vars
        self.fluid_vars += self.coupling_fluid_vars

        if parameters["pressure_grad_recovery_type"].GetInt() > 0:
            self.fluid_vars += [Fluid.RECOVERED_PRESSURE_GRADIENT]

        if (parameters["gradient_calculation_type"].GetInt() > 1
            or parameters["pressure_grad_recovery_type"].GetInt() > 1
            or parameters["material_acceleration_calculation_type"].GetInt() == 7
            or parameters["laplacian_calculation_type"].GetInt() > 1
            or parameters["coupling"]["backward_coupling"]["fluid_fraction_grad_type"] == 3):
            self.fluid_vars += [Fluid.NODAL_WEIGHTS]

        if parameters["material_acceleration_calculation_type"].GetInt():
            self.fluid_vars += [Kratos.MATERIAL_ACCELERATION]
            self.fluid_vars += [Kratos.VELOCITY_COMPONENT_GRADIENT]

            # if (parameters["material_acceleration_calculation_type"].GetInt() == 5
            #     or parameters["material_acceleration_calculation_type"].GetInt() == 6):
            if parameters["store_full_gradient_option"].GetBool():
                self.fluid_vars += [Kratos.VELOCITY_X_GRADIENT]
                self.fluid_vars += [Kratos.VELOCITY_Y_GRADIENT]
                self.fluid_vars += [Kratos.VELOCITY_Z_GRADIENT]

        if (parameters["vorticity_calculation_type"].GetInt() > 0
            or PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')):
            self.fluid_vars += [Kratos.VORTICITY]

        if parameters["laplacian_calculation_type"].GetInt():
            self.fluid_vars += [Kratos.VELOCITY_LAPLACIAN]

        if PT.RecursiveFindTrueBoolInParameters(parameters["properties"], 'do_apply_faxen_corrections'):
            self.fluid_vars += [Kratos.VELOCITY_LAPLACIAN_RATE]

        if parameters["coupling"]["backward_coupling"]["calculate_diffusivity_option"].GetBool():
            self.fluid_vars += [Kratos.CONDUCTIVITY]

        # dem variables
        self.dem_vars = []
        self.dem_vars += self.dem_printing_vars
        self.dem_vars += self.coupling_dem_vars
        self.dem_vars += [Kratos.BUOYANCY]
        self.dem_vars += [Kratos.VELOCITY_OLD]
        self.dem_vars += [Kratos.SLIP_VELOCITY]
        self.dem_vars += [DEM.DEM_STRESS_TENSOR]
        self.dem_vars += [DEM.GROUP_ID]

        if self.do_include_history_force:
            self.dem_vars += [Kratos.BASSET_FORCE]

        if parameters["frame_of_reference"]["frame_type"].GetInt() and self.do_include_history_force > 0:
            self.dem_vars += [SDEM.DISPLACEMENT_OLD]
            self.dem_vars += [Kratos.VELOCITY_OLD_OLD]

        if (parameters["custom_dem"]["translational_integration_scheme"].GetString()
            in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}
            or self.do_include_history_force):
            self.dem_vars += [Kratos.VELOCITY_OLD]
            self.dem_vars += [Kratos.ADDITIONAL_FORCE_OLD]
            self.dem_vars += [Kratos.AUX_VEL]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [Kratos.DRAG_FORCE]

        self.dem_vars += [Kratos.PARTICLE_SPHERICITY] # TODO: add only when needed

        if (PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')
            and parameters["add_each_hydro_force_option"].GetBool()):
            self.dem_vars += [Kratos.LIFT_FORCE]

        # if parameters["add_each_hydro_force_option"].GetBool() and parameters["properties"][0]["hydrodynamic_law_parameters"]["virtual_mass_force_parameters"]["name"].GetString() != "default":
        #     self.dem_vars += [Kratos.VIRTUAL_MASS_FORCE]

        # if parameters["add_each_hydro_force_option"].GetBool() and parameters["properties"][0]["hydrodynamic_law_parameters"]["undisturbed_force_parameters"]["name"].GetString() != "default":
        #     self.dem_vars += [SDEM.UNDISTURBED_FLOW_FORCE]

        # clusters variables
        self.clusters_vars = []

        # rigid faces variables
        self.rigid_faces_vars = [Kratos.VELOCITY,
                                 Kratos.ANGULAR_VELOCITY,
                                 Kratos.DISPLACEMENT,
                                 DEM.DELTA_DISPLACEMENT,
                                 Kratos.DELTA_ROTATION,
                                 DEM.CONTACT_FORCES,
                                 DEM.DEM_PRESSURE,
                                 DEM.ELASTIC_FORCES,
                                 Kratos.PRESSURE,
                                 DEM.TANGENTIAL_ELASTIC_FORCES,
                                 DEM.SHEAR_STRESS,
                                 Kratos.NODAL_AREA,
                                 Kratos.VELOCITY_OLD]

        if parameters["custom_fluid"]["embedded_option"].GetBool():
            self.rigid_faces_vars += [Kratos.FORCE]
            self.rigid_faces_vars += [Kratos.POSITIVE_FACE_PRESSURE]
            self.rigid_faces_vars += [Kratos.NEGATIVE_FACE_PRESSURE]

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

        if parameters["custom_fluid"]["embedded_option"].GetBool():
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
            self.fluid_printing_vars += [GetGlobalVariableByName(variable)]

        for variable in self.dem_nodal_results:
            self.dem_printing_vars += [GetGlobalVariableByName(variable)]

        for variable in self.clusters_nodal_results:
            self.clusters_printing_vars += [GetGlobalVariableByName(variable)]

        for variable in self.rigid_faces_nodal_results:
            self.rigid_faces_printing_vars += [GetGlobalVariableByName(variable)]

        for variable in self.mixed_nodal_results:
            self.dem_printing_vars += [GetGlobalVariableByName(variable)]
            self.fluid_printing_vars += [GetGlobalVariableByName(variable)]

        for var in self.mixed_nodal_results:

            if var in self.nodal_results:
                self.nodal_results.remove(var)

        VariablesManager.EliminateRepeatedValuesFromList(self.nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.dem_nodal_results)
        VariablesManager.EliminateRepeatedValuesFromList(self.mixed_nodal_results)

    def ConstructListsOfVariablesForCoupling(self, parameters):

        # fluid coupling variables
        self.coupling_fluid_vars = []
        self.coupling_fluid_vars += [Kratos.MATERIAL_ACCELERATION]
        self.coupling_fluid_vars += [Fluid.MASS_SOURCE]
        #self.coupling_fluid_vars += [Kratos.EXACT_VELOCITY]
        # self.coupling_fluid_vars += [SDEM.VECTORIAL_ERROR]
        self.coupling_fluid_vars += [SDEM.ERROR_X]
        self.coupling_fluid_vars += [SDEM.ERROR_Y]
        self.coupling_fluid_vars += [SDEM.ERROR_Z]
        self.coupling_fluid_vars += [SDEM.ERROR_P]
        self.coupling_fluid_vars += [SDEM.SCALAR_ERROR]
        self.coupling_fluid_vars += [SDEM.EXACT_PRESSURE]

        self.coupling_fluid_vars += [Kratos.KratosGlobals.GetVariable(parameters["body_force_per_unit_mass_variable_name"].GetString() )]

        # Nodal density
        if 'EXACT_MATERIAL_ACCELERATION' in self.dem_nodal_results:
            self.coupling_fluid_vars += [SDEM.EXACT_MATERIAL_ACCELERATION]

        if parameters["custom_fluid"]["fluid_model_type"].GetInt() == 0:
            self.coupling_fluid_vars += [SDEM.AVERAGED_FLUID_VELOCITY]

        if (parameters["custom_fluid"]["fluid_model_type"].GetInt() == 0
            or parameters["coupling"]["coupling_level_type"].GetInt() >= 1):
            self.coupling_fluid_vars += [Kratos.FLUID_FRACTION]
            self.coupling_fluid_vars += [Kratos.FLUID_FRACTION_OLD]
            self.coupling_fluid_vars += [SDEM.FLUID_FRACTION_OLD_2]

            if 'DISPERSE_FRACTION' in self.nodal_results:
                self.coupling_fluid_vars += [Kratos.DISPERSE_FRACTION]

            self.coupling_fluid_vars += [Kratos.AVERAGED_PARTICLE_VELOCITY]
            if parameters["coupling"]["backward_coupling"]["filter_velocity_option"].GetBool():
                self.coupling_fluid_vars += [Kratos.PARTICLE_VEL_FILTERED]
                self.coupling_fluid_vars += [Kratos.TIME_AVERAGED_ARRAY_3]
                self.coupling_fluid_vars += [Kratos.PHASE_FRACTION]

        if parameters["custom_fluid"]["fluid_model_type"].GetInt() >= 1:
            self.coupling_fluid_vars += [Kratos.FLUID_FRACTION_GRADIENT]
            self.coupling_fluid_vars += [Kratos.FLUID_FRACTION_RATE]

        if parameters["coupling"]["coupling_level_type"].GetInt() >= 1:
            self.coupling_fluid_vars += [Kratos.HYDRODYNAMIC_REACTION]

        if parameters["coupling"]["coupling_level_type"].GetInt() >= 1 and parameters["coupling"]["forward_coupling"]["time_averaging_type"].GetInt() > 0:
            self.coupling_fluid_vars += [Kratos.MEAN_HYDRODYNAMIC_REACTION]

        if parameters["non_newtonian_fluid"]["non_newtonian_option"].GetBool():
            self.coupling_fluid_vars += [Kratos.POWER_LAW_N]
            self.coupling_fluid_vars += [Kratos.POWER_LAW_K]
            self.coupling_fluid_vars += [Kratos.YIELD_STRESS]
            # self.coupling_fluid_vars += [Kratos.GEL_STRENGTH] # TODO: make specific option for this

        if parameters["coupling"]["backward_coupling"]["viscosity_modification_type"].GetInt():
            self.coupling_fluid_vars += [Kratos.VISCOSITY]

        if parameters["custom_fluid"]["embedded_option"].GetBool():
            self.coupling_fluid_vars += [Kratos.DISTANCE]

        # dem coupling variables
        self.coupling_dem_vars = []

        # Nodal density
        if 'NODAL_DENSITY_PROJECTED' in self.dem_nodal_results:
            self.coupling_dem_vars += [SDEM.NODAL_DENSITY_PROJECTED]
        if 'TRUNC_SOLUTION_PROJECTED' in self.dem_nodal_results:
            self.coupling_dem_vars += [SDEM.TRUNC_SOLUTION_PROJECTED]
        if 'NODAL_DENSITY' in self.dem_nodal_results:
            self.coupling_fluid_vars += [SDEM.NODAL_DENSITY]
        if 'TRUNC_SOLUTION_SDEM' in self.dem_nodal_results:
            self.coupling_fluid_vars += [SDEM.TRUNC_SOLUTION_SDEM]
        if 'EXACT_MATERIAL_ACCELERATION' in self.dem_nodal_results:
            self.coupling_dem_vars += [SDEM.EXACT_MATERIAL_ACCELERATION]

        if parameters["coupling"]["coupling_level_type"].GetInt() > 0:
            self.coupling_dem_vars += [Kratos.FLUID_VEL_PROJECTED]
            self.coupling_dem_vars += [Kratos.PRESSURE_GRAD_PROJECTED]
            self.coupling_dem_vars += [Kratos.FLUID_ACCEL_PROJECTED]
            self.coupling_dem_vars += [Kratos.FLUID_DENSITY_PROJECTED]
            self.coupling_dem_vars += [Kratos.DRAG_COEFFICIENT]
            self.coupling_dem_vars += [SDEM.HYDRODYNAMIC_REACTION_PROJECTED]
            self.coupling_dem_vars += [Kratos.FLUID_VISCOSITY_PROJECTED]
            self.coupling_dem_vars += [Kratos.HYDRODYNAMIC_FORCE]
            self.coupling_dem_vars += [Kratos.HYDRODYNAMIC_MOMENT]
            self.coupling_dem_vars += [Kratos.MATERIAL_FLUID_ACCEL_PROJECTED]
            self.coupling_dem_vars += [Kratos.FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED]
            self.coupling_dem_vars += [Kratos.ADDITIONAL_FORCE] # Here for safety for the moment

            if PT.RecursiveFindTrueBoolInParameters(parameters["properties"], 'do_apply_faxen_corrections'):
                self.coupling_dem_vars += [Kratos.FLUID_VEL_LAPL_PROJECTED]
                self.coupling_dem_vars += [Kratos.FLUID_VEL_LAPL_RATE_PROJECTED]

            if self.do_include_history_force:
                self.coupling_dem_vars += [Kratos.FLUID_VEL_PROJECTED_RATE]

        if parameters["coupling"]["coupling_level_type"].GetInt() >= 1 or parameters["custom_fluid"]["fluid_model_type"].GetInt() == 0:
            self.coupling_dem_vars += [Kratos.FLUID_FRACTION_PROJECTED]

        if (parameters["coupling"]["coupling_level_type"].GetInt() >= 1
            and 'FLUID_FRACTION_GRADIENT_PROJECTED' in self.dem_printing_vars):
            self.coupling_dem_vars += [Kratos.FLUID_FRACTION_GRADIENT_PROJECTED]

        if parameters["non_newtonian_fluid"]["non_newtonian_option"].GetBool():
            self.coupling_dem_vars += [Kratos.POWER_LAW_N]
            self.coupling_dem_vars += [Kratos.POWER_LAW_K]
            self.coupling_dem_vars += [Kratos.YIELD_STRESS]

        if PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters', condition=lambda value: not (value['name'].GetString()=='default')):
            self.coupling_dem_vars += [Kratos.FLUID_VORTICITY_PROJECTED]
            self.coupling_dem_vars += [Kratos.SHEAR_RATE_PROJECTED]

        if parameters["custom_fluid"]["embedded_option"].GetBool():
            self.coupling_dem_vars += [Kratos.DISTANCE]

        if 'REYNOLDS_NUMBER' in self.dem_nodal_results:
            self.coupling_dem_vars += [Kratos.REYNOLDS_NUMBER]

        if self.do_backward_coupling:
            self.coupling_dem_vars += [SDEM.GENTLE_INITIATION_COUPLING_COEFFICIENT]

            if parameters["coupling"]["backward_coupling"]["apply_time_filter_to_fluid_fraction_option"].GetBool():
                self.time_filtered_vars += [Kratos.FLUID_FRACTION]

        if parameters["coupling"]["backward_coupling"]["filter_velocity_option"].GetBool():
            self.coupling_fluid_vars += [Kratos.AVERAGED_PARTICLE_VELOCITY]
            self.time_filtered_vars += [Kratos.PARTICLE_VEL_FILTERED]

        if self.time_filtered_vars:
            self.coupling_fluid_vars += [Kratos.TIME_AVERAGED_DOUBLE]
            self.coupling_fluid_vars += [Kratos.TIME_AVERAGED_ARRAY_3]
            self.coupling_fluid_vars += [SDEM.TIME_AVERAGED_BODY_FORCE]
            self.time_filtered_vars  += [Kratos.KratosGlobals.GetVariable(parameters["body_force_per_unit_mass_variable_name"].GetString())]

    def ChangeListOfFluidNodalResultsToPrint(self, parameters):
        fluid_list = self.project_parameters["fluid_nodal_results"]
        self.nodal_results.extend(key for key in fluid_list.keys() if fluid_list[key].GetBool())

        if parameters["store_full_gradient_option"].GetBool() and 'VELOCITY_GRADIENT' in self.nodal_results:
            self.nodal_results += ["VELOCITY_X_GRADIENT"]
            self.nodal_results += ["VELOCITY_Y_GRADIENT"]
            self.nodal_results += ["VELOCITY_Z_GRADIENT"]

        if parameters["custom_fluid"]["fluid_model_type"].GetInt() == 0 and 'AVERAGED_FLUID_VELOCITY' in self.nodal_results:
            self.nodal_results += ["AVERAGED_FLUID_VELOCITY"]

        if parameters["custom_fluid"]["fluid_model_type"].GetInt() == 1 and 'FLUID_FRACTION_GRADIENT' in self.nodal_results:
            self.nodal_results += ["FLUID_FRACTION_GRADIENT"]
