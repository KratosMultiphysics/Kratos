import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def _DotProduct(A,B):
    return sum(i[0]*i[1] for i in zip(A, B))

def _CrossProduct(A, B):
    C = KratosMultiphysics.Vector(3)
    C[0] = A[1]*B[2]-A[2]*B[1]
    C[1] = A[2]*B[0]-A[0]*B[2]
    C[2] = A[0]*B[1]-A[1]*B[0]
    return C

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "please specify the model part that contains the surface nodes",
            "far_field_model_part_name": "please specify the model part that contains the surface nodes",
            "moment_reference_point" : [0.0,0.0,0.0]
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        far_field_model_part_name = settings["far_field_model_part_name"].GetString()
        if far_field_model_part_name == "":
            err_msg = "Empty model_part_name in ComputeLiftProcess\n"
            err_msg += "Please specify the model part that contains the far field surface nodes"
            raise Exception(err_msg)
        self.far_field_model_part = Model[far_field_model_part_name]
        self.fluid_model_part = self.body_model_part.GetRootModelPart()
        self.reference_area =  self.fluid_model_part.ProcessInfo.GetValue(CPFApp.REFERENCE_CHORD)
        self.moment_reference_point = settings["moment_reference_point"].GetVector()

        if not self.reference_area > 0.0:
            raise Exception('The reference area should be larger than 0.')

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','COMPUTE LIFT')

        self.__CalculateWakeTangentAndNormalDirections()
        self.__ComputeLiftFromPressure()
        self.__ComputeLiftFromFarField()
        if(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.__ComputeMomentFromPressure()
            self.__ComputeLiftFromJumpCondition()

    def __CalculateWakeTangentAndNormalDirections(self):
        self.wake_direction = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)
        dnorm = self.wake_direction.norm_2()
        self.wake_direction /= dnorm

        self.wake_normal = KratosMultiphysics.Vector(3)
        if(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.wake_normal[0] = -self.wake_direction[1]
            self.wake_normal[1] = self.wake_direction[0]
            self.wake_normal[2] = 0.0
        elif(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            # TODO: Read wake normal from wake process
            self.wake_normal[0] = 0.0
            self.wake_normal[1] = 0.0
            self.wake_normal[2] = 1.0

        self.span_direction = KratosMultiphysics.Vector(3)
        self.span_direction = _CrossProduct(self.wake_normal, self.wake_direction)

    def __ComputeLiftFromPressure(self):
        force_coefficient = KratosMultiphysics.Vector(3)
        self.moment_coefficient = KratosMultiphysics.Vector(3)

        for cond in self.body_model_part.Conditions:
            surface_normal = cond.GetGeometry().Normal()
            pressure_coefficient = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)

            # Computing forces
            force_coefficient += surface_normal*pressure_coefficient

        force_coefficient /= self.reference_area

        self.lift_coefficient = _DotProduct(force_coefficient,self.wake_normal)
        self.drag_coefficient = _DotProduct(force_coefficient,self.wake_direction)
        self.lateral_force_coefficient = _DotProduct(force_coefficient,self.span_direction)

        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cl = ', self.lift_coefficient)
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cd = ', self.drag_coefficient)
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cq = ', self.lateral_force_coefficient)

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT, self.lift_coefficient)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.DRAG_COEFFICIENT, self.drag_coefficient)

    def __ComputeMomentFromPressure(self):
        self.moment_coefficient = KratosMultiphysics.Vector(3)

        for cond in self.body_model_part.Conditions:
            surface_normal = cond.GetGeometry().Normal()
            pressure_coefficient = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)

            # Computing moment
            mid_point = cond.GetGeometry().Center()
            lever = mid_point-self.moment_reference_point
            self.moment_coefficient += _CrossProduct(lever, surface_normal*(-pressure_coefficient))

        self.moment_coefficient /= self.reference_area

        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cm = ', self.moment_coefficient[2])
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.MOMENT_COEFFICIENT, self.moment_coefficient[2])

    def __ComputeLiftFromJumpCondition(self):
        # Find the Trailing Edge node
        for node in self.body_model_part.Nodes:
            if node.GetValue(CPFApp.TRAILING_EDGE):
                te = node
                break

        free_stream_velocity = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)
        u_inf = free_stream_velocity.norm_2()

        node_velocity_potential_te = te.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
        node_auxiliary_velocity_potential_te = te.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
        if(te.GetValue(CPFApp.WAKE_DISTANCE) > 0.0):
            potential_jump_phi_minus_psi_te = node_velocity_potential_te - node_auxiliary_velocity_potential_te
        else:
            potential_jump_phi_minus_psi_te = node_auxiliary_velocity_potential_te - node_velocity_potential_te
        self.lift_coefficient_jump = 2*potential_jump_phi_minus_psi_te/ ( u_inf * self.reference_area )

        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cl = ' , self.lift_coefficient_jump, ' = ( 2 * DPhi ) / ( U_inf * c )')
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT_JUMP, self.lift_coefficient_jump)

    def __ComputeLiftFromFarField(self):

        force_coefficient_pres = KratosMultiphysics.Vector(3)
        force_coefficient_vel = KratosMultiphysics.Vector(3)
        free_stream_velocity = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)

        for cond in self.far_field_model_part.Conditions:
            surface_normal = cond.GetGeometry().Normal()

            # Computing contribution due to pressure
            pressure_coefficient = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            force_coefficient_pres -= surface_normal*pressure_coefficient

            # Computing contribution due to convection
            velocity = cond.GetValue(KratosMultiphysics.VELOCITY)
            velocity_projection = _DotProduct(velocity, surface_normal)
            disturbance = velocity - free_stream_velocity
            density = cond.GetValue(KratosMultiphysics.DENSITY)
            force_coefficient_vel -= velocity_projection * disturbance * density

        # Normalizing with reference area
        force_coefficient_pres /= self.reference_area

        # Normalizing with dynamic pressure and reference area
        free_stream_density = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_DENSITY)
        free_stream_velocity_norm2 = _DotProduct(free_stream_velocity,free_stream_velocity)
        dynamic_pressure = 0.5 * free_stream_velocity_norm2 * free_stream_density
        force_coefficient_vel /= dynamic_pressure * self.reference_area

        force_coefficient = force_coefficient_pres + force_coefficient_vel
        self.lift_coefficient_far_field = _DotProduct(force_coefficient,self.wake_normal)
        self.drag_coefficient_far_field = _DotProduct(force_coefficient,self.wake_direction)

        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cl = ', self.lift_coefficient_far_field, 'Far field')
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess',' Cd = ', self.drag_coefficient_far_field, 'Far field')

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT_FAR_FIELD, self.lift_coefficient_far_field)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.DRAG_COEFFICIENT_FAR_FIELD, self.drag_coefficient_far_field)
