import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def _DotProduct(A,B):
    return sum(i[0]*i[1] for i in zip(A, B))

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeFarFieldLiftProcess(Model, settings["Parameters"])

class ComputeFarFieldLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "please specify the model part that contains the far field nodes"
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.far_field_model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.far_field_model_part.GetRootModelPart()
        self.reference_area =  self.fluid_model_part.ProcessInfo.GetValue(CPFApp.REFERENCE_CHORD)

        if not self.reference_area > 0.0:
            raise Exception('The reference area should be larger than 0.')

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeFarFieldLiftProcess','COMPUTE LIFT')

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

        # Normalizing with static pressure and reference area
        free_stream_density = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_DENSITY)
        free_stream_velocity_norm2 = _DotProduct(free_stream_velocity,free_stream_velocity)
        static_pressure = 0.5 * free_stream_velocity_norm2 * free_stream_density
        force_coefficient_vel /= static_pressure * self.reference_area

        self.__CalculateWakeTangentAndNormalDirections()

        force_coefficient = force_coefficient_pres + force_coefficient_vel
        self.lift_coefficient = _DotProduct(force_coefficient,self.wake_normal)
        self.drag_coefficient = _DotProduct(force_coefficient,self.wake_direction)

        KratosMultiphysics.Logger.PrintInfo('ComputeFarFieldLiftProcess',' Cl = ', self.lift_coefficient)
        KratosMultiphysics.Logger.PrintInfo('ComputeFarFieldLiftProcess',' Cd = ', self.drag_coefficient)

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT_FAR_FIELD, self.lift_coefficient)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.DRAG_COEFFICIENT_FAR_FIELD, self.drag_coefficient)

    def __CalculateWakeTangentAndNormalDirections(self):
        self.wake_direction = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)
        if(self.wake_direction.Size() != 3):
            raise Exception('The wake direction should be a vector with 3 components!')

        dnorm = self.wake_direction.norm_2()
        self.wake_direction /= dnorm

        self.wake_normal = KratosMultiphysics.Vector(3)
        self.wake_normal[0] = -self.wake_direction[1]
        self.wake_normal[1] = self.wake_direction[0]
        self.wake_normal[2] = 0.0
