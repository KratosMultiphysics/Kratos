import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldAndWakeProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFarFieldAndWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        if not settings.Has("angle_of_attack_units"):
            KratosMultiphysics.Logger.PrintWarning("ApplyFarFieldAndWakeProcess", "'angle of attack_units' is not provided. Using 'radians' as default angle of attack unit.")

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"",
                "angle_of_attack_units": "radians",
                "angle_of_attack": 0.0,
                "mach_infinity": 0.02941176471,
                "free_stream_density"  : 1.0,
                "speed_of_sound": 340,
                "heat_capacity_ratio": 1.4,
                "inlet_potential": 1.0,
                "mach_number_limit": 0.94,
                "mach_number_squared_limit": -1,
                "critical_mach": 0.99,
                "upwind_factor_constant": 1.0,
                "initialize_flow_field": true,
                "perturbation_field": false,
                "define_wake": false,
                "check_wake_condition_tolerance" : 1e-1,
                "check_wake_condition_echo_level" : 0,
                "compute_wake_at_each_step"      : false,
                "wake_type"  : "",
                "wake_parameters": {
                }
            }  """ )
        settings.ValidateAndAssignDefaults(default_parameters)


        self.far_field_model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.far_field_model_part.GetRootModelPart()

        self.angle_of_attack_units = settings["angle_of_attack_units"].GetString()
        self.angle_of_attack = settings["angle_of_attack"].GetDouble()
        self.free_stream_mach = settings["mach_infinity"].GetDouble()
        self.density_inf = settings["free_stream_density"].GetDouble()
        self.free_stream_speed_of_sound = settings["speed_of_sound"].GetDouble()
        self.heat_capacity_ratio = settings["heat_capacity_ratio"].GetDouble()
        self.inlet_potential_0 = settings["inlet_potential"].GetDouble()
        self.mach_number_limit = settings["mach_number_limit"].GetDouble()
        self.mach_number_squared_limit = settings["mach_number_squared_limit"].GetDouble()
        self.critical_mach = settings["critical_mach"].GetDouble()
        self.upwind_factor_constant = settings["upwind_factor_constant"].GetDouble()
        self.initialize_flow_field = settings["initialize_flow_field"].GetBool()
        self.perturbation_field = settings["perturbation_field"].GetBool()
        self.define_wake = settings["define_wake"].GetBool()
        self.check_wake_condition_tolerance = settings["check_wake_condition_tolerance"].GetDouble()
        self.compute_wake_at_each_step = settings["compute_wake_at_each_step"].GetBool()
        self.check_wake_condition_echo_level = settings["check_wake_condition_echo_level"].GetInt()

        if(self.perturbation_field):
            self.initialize_flow_field = False

        # Computing free stream velocity
        self.u_inf = self.free_stream_mach * self.free_stream_speed_of_sound
        self.free_stream_velocity = KratosMultiphysics.Vector(3)

        if self.angle_of_attack_units == "radians":
            KratosMultiphysics.Logger.PrintWarning("ApplyFarFieldAndWakeProcess", " Using 'radians' as angle of attack unit. This will be deprecated soon.")
        elif self.angle_of_attack_units == "degrees":
            self.angle_of_attack = self.angle_of_attack*math.pi/180

        self.domain_size = self.fluid_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        if self.domain_size == 2:
            # By convention 2D airfoils are in the xy plane
            self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
            self.free_stream_velocity[1] = round(self.u_inf*math.sin(self.angle_of_attack),8)
            self.free_stream_velocity[2] = 0.0
        else: # self.domain_size == 3
            # By convention 3D wings and aircrafts:
            # y axis along the span
            # z axis pointing upwards
            # TODO: Add sideslip angle beta
            self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
            self.free_stream_velocity[1] = 0.0
            self.free_stream_velocity[2] = round(self.u_inf*math.sin(self.angle_of_attack),8)

        self.free_stream_velocity_direction = self.free_stream_velocity / self.u_inf

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_MACH,self.free_stream_mach)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_VELOCITY,self.free_stream_velocity)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_VELOCITY_DIRECTION,self.free_stream_velocity_direction)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_DENSITY,self.density_inf)
        self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.SOUND_VELOCITY,self.free_stream_speed_of_sound)
        self.fluid_model_part.ProcessInfo.SetValue(KratosCFD.HEAT_CAPACITY_RATIO,self.heat_capacity_ratio)

        if self.mach_number_squared_limit > 0.0:
            self.fluid_model_part.ProcessInfo.SetValue(CPFApp.MACH_LIMIT,math.sqrt(self.mach_number_squared_limit))
            warn_msg = 'Both mach_number_squared_limit and mach_number_limit are defined. Using mach_number_squared_limit = ' + str(self.mach_number_squared_limit)
            KratosMultiphysics.Logger.PrintWarning('ApplyFarFieldAndWakeProcess', warn_msg)
        else:
            self.fluid_model_part.ProcessInfo.SetValue(CPFApp.MACH_LIMIT,self.mach_number_limit)

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.CRITICAL_MACH,self.critical_mach)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.UPWIND_FACTOR_CONSTANT,self.upwind_factor_constant)

        # Setting wake operation
        if self.define_wake:
            registry_entry = settings["wake_type"].GetString()
            operation = KratosMultiphysics.Registry[f"{registry_entry}.Prototype"]
            self.wake_operation = operation.Create(Model, settings["wake_parameters"])

    def ExecuteInitialize(self):
        far_field_process=CPFApp.ApplyFarFieldProcess(self.far_field_model_part, self.inlet_potential_0, self.initialize_flow_field, self.perturbation_field)
        far_field_process.Execute()

        if self.define_wake:
            self.wake_operation.Execute()

    def ExecuteInitializeSolutionStep(self):
        # TODO REMOVE FROM HERE, Used only for optimization
        if self.define_wake:
            if  self.compute_wake_at_each_step and self.fluid_model_part.ProcessInfo[KratosMultiphysics.STEP] > 1:
                KratosMultiphysics.Logger.PrintWarning("ApplyFarFieldAndWakeProcess", " Using define_wake_operation at each step. This will be deprecated soon.")
                self.wake_operation.Execute()

    def ExecuteFinalizeSolutionStep(self):
        if self.define_wake:
            if not self.fluid_model_part.HasSubModelPart("wake_elements_model_part"):
                raise Exception("Fluid model part does not have a wake_elements_model_part")
            else: self.wake_elements_model_part = self.fluid_model_part.GetSubModelPart("wake_elements_model_part")
            
            if self.domain_size == 2:
                CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled2D(
                    self.wake_elements_model_part, self.check_wake_condition_tolerance, self.check_wake_condition_echo_level)
                CPFApp.PotentialFlowUtilities.ComputePotentialJump2D(self.wake_elements_model_part)
            else: # self.domain_size == 3
                CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(
                    self.wake_elements_model_part, self.check_wake_condition_tolerance, self.check_wake_condition_echo_level)
                CPFApp.PotentialFlowUtilities.ComputePotentialJump3D(self.wake_elements_model_part)
