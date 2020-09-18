import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignGravityToParticleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignGravityToParticleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "MP_VOLUME_ACCELERATION",
                "modulus"              : 1.0,
                "constrained"          : true,
                "direction"            : [0.0, 0.0, 0.0],
                "local_axes"           : {}
            }
            """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("modulus")):
            if(settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")

        if(settings.Has("direction")):
            if(settings["direction"].IsString()):
                default_settings["direction"].SetString("Automatic")

        # Detect if variable_name is MP_VOLUME_ACCELERATION
        if(settings.Has("variable_name")):
            if(settings["variable_name"].GetString() != "MP_VOLUME_ACCELERATION"):
                KratosMultiphysics.Logger.PrintInfo("Warning in apply gravity to particle", "Error in determining variable_name")
                raise Exception('The assign_gravity_to_particle_process only accepts \"MP_VOLUME_ACCELERATION\" as variable_name.')

        settings.ValidateAndAssignDefaults(default_settings)

        # Default settings
        self.model_part = Model[settings["model_part_name"].GetString()]
        
        self.gravity_direction = settings["direction"].GetVector()
        #self.gravity_acceleration = self.modulus * self.gravity_direction
        self.value_is_numeric = False
        if settings["modulus"].IsNumber():
            self.value_is_numeric = True
            self.modulus = settings["modulus"].GetDouble()
            self.gravity_acceleration = self.gravity_direction * self.modulus
        else:
            self.function_string = settings["modulus"].GetString()
            self.aux_function = KratosMultiphysics.PythonGenericFunctionUtility(self.function_string, settings["local_axes"])


    def ExecuteBeforeSolutionLoop(self):
        # Assign gravity to MP after solver.Initialize() - only apply once at the beginning!
        for element in self.model_part.Elements:
            if self.value_is_numeric:
                self.gravity_acceleration = self.gravity_direction * self.modulus
                element.SetValuesOnIntegrationPoints(KratosParticle.MP_VOLUME_ACCELERATION,[self.gravity_acceleration],self.model_part.ProcessInfo)
                element.SetValuesOnIntegrationPoints(KratosParticle.MP_ACCELERATION,[self.gravity_acceleration],self.model_part.ProcessInfo)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        
        for element in self.model_part.Elements:
            step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            time_step = current_time / step
            
           
            if not self.value_is_numeric:
                if self.aux_function.DependsOnSpace() == False: #depends on time only
                    delta_gravity_acceleration = self.aux_function.CallFunction(0.0,0.0,0.0,time_step,0.0,0.0,0.0) * self.gravity_direction

                volume_gravity_acceleration = element.CalculateOnIntegrationPoints(KratosParticle.MP_VOLUME_ACCELERATION,self.model_part.ProcessInfo)[0]
                gravity_acceleration = element.CalculateOnIntegrationPoints(KratosParticle.MP_ACCELERATION,self.model_part.ProcessInfo)[0]

                volume_gravity_acceleration += delta_gravity_acceleration
                gravity_acceleration += delta_gravity_acceleration
                        
                
                element.SetValuesOnIntegrationPoints(KratosParticle.MP_VOLUME_ACCELERATION,[volume_gravity_acceleration],self.model_part.ProcessInfo)
                element.SetValuesOnIntegrationPoints(KratosParticle.MP_ACCELERATION,[gravity_acceleration],self.model_part.ProcessInfo)
                
            