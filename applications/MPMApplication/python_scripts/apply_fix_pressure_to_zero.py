
import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFixPressureToZero(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFixPressureToZero(KratosMultiphysics.Process):
        def __init__(self, Model, settings ):
            KratosMultiphysics.Process.__init__(self)

            default_settings = KratosMultiphysics.Parameters("""
                {
                    "mesh_id"              : 0,
                    "model_part_name"      : "please_specify_model_part_name",
                    "variable_name"        : "PRESSURE",
                    "modulus"              : 1.0,
                    "local_axes"           : {}
                }
                """)

            # Trick: allow "modulus" and "component" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)

            if settings.Has("modulus"):
                if settings["modulus"].IsString():
                    default_settings["modulus"].SetString("0.0")


            self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
    

            settings.ValidateAndAssignDefaults(default_settings)


            self.vector_direction = KratosMultiphysics.Vector(3)
            self.aux_function_direction = ["0.0","0.0","0.0"]
            self.value_direction_is_numeric = [False, False, False]
       
            self.value_is_numeric = False
            self.value = KratosMultiphysics.Vector(3)
            self.acceleration = KratosMultiphysics.Vector(3)
            self.aux_function = "0.0"

            if settings["modulus"].IsNumber():
                self.value_is_numeric = True
                self.modulus = settings["modulus"].GetDouble()
                self.value = self.vector_direction * self.modulus
            else:
                self.function_string = settings["modulus"].GetString()
                self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])

            # Default settings
            self.model_part = Model[settings["model_part_name"].GetString()]


        def ExecuteInitializeSolutionStep(self):                
                for mp in self.model_part.elements:
                    mp_coord = mp.CalculateOnIntegrationPoints(KratosMPM.MP_COORD,self.model_part.ProcessInfo)[0]
                    if ((mp_coord[0]<-0.49) and (mp_coord[1]<0.0001)):
                        mp.Fix(KratosMultiphysics.PRESSURE)
                        mp.SetSolutionStepValue(KratosMultiphysics.MP_PRESSURE, 0, 0)  