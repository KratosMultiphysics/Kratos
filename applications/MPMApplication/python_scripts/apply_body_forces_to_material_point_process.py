import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignBodyForcesToMaterialPointProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignBodyForcesToMaterialPointProcess(KratosMultiphysics.Process):
        def __init__(self, Model, settings ):
            KratosMultiphysics.Process.__init__(self)

            default_settings = KratosMultiphysics.Parameters("""
                {
                    "mesh_id"              : 0,
                    "model_part_name"      : "please_specify_model_part_name",
                    "variable_name"        : "MP_VOLUME_ACCELERATION",
                    "modulus"              : 1.0,
                    "constrained"          : true,
                    "component"            : [0.0, 0.0, 0.0],
                    "local_axes"           : {}
                }
                """)

            # Trick: allow "modulus" and "component" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)

            if settings.Has("modulus"):
                if settings["modulus"].IsString():
                    default_settings["modulus"].SetString("0.0")

            if settings.Has("component"):
                if settings["component"].IsString():
                    default_settings["component"].SetString("Automatic")

            self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
            # Detect if variable_name is MP_VOLUME_ACCELERATION
            #if(settings.Has("variable_name")):
            #    if(settings["variable_name"].GetString() != "MP_VOLUME_ACCELERATION"):
            #        KratosMultiphysics.Logger.PrintInfo("Warning in apply gravity to material point", "Error in determining variable_name")
            #        raise Exception('The assign_gravity_to_material_point_process only accepts \"MP_VOLUME_ACCELERATION\" as variable_name.')

            settings.ValidateAndAssignDefaults(default_settings)


            self.vector_direction = KratosMultiphysics.Vector(3)
            self.aux_function_direction = ["0.0","0.0","0.0"]
            self.value_direction_is_numeric = [False, False, False]
            for i in range(3):
                if settings["component"][i].IsNumber():
                    self.value_direction_is_numeric[i] = True
                    self.vector_direction[i] = settings["component"][i].GetDouble()
                else:
                    self.function_string_direction = settings["component"][i].GetString()
                    self.aux_function_direction[i] = KratosMultiphysics.GenericFunctionUtility(self.function_string_direction, settings["local_axes"])


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
            #model_part = self.model[self.mpm_material_model_part_name]
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

            # the model part is identified here, AFTER it has been transferred to the MPM_material part!
            if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                 #for node in self.model_part.Nodes:
                 for mp in self.model_part.Elements:
                    mp_coord = mp.CalculateOnIntegrationPoints(KratosMPM.MP_COORD,self.model_part.ProcessInfo)[0]
                    

                    for i in range(3):
                        if not self.value_direction_is_numeric[i]:
                            if self.aux_function_direction[i].DependsOnSpace() == False: #depends on time only
                                self.vector_direction[i] = self.aux_function_direction[i].CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                                
                            #else: #most general case - space varying function (possibly also time varying)
                                #self.vector_direction[i]= self.aux_function_direction[i].CallFunction(mp_coord[0],mp_coord[1],mp_coord[2],current_time,0.0,0.0,0.0)


                        if not self.value_is_numeric:
                            if self.aux_function.DependsOnSpace() == False: #depends on time only
                                self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0) * self.vector_direction
                                 # - space varying function (possibly also time varying)
                                #self.value = self.aux_function.CallFunction(mp_coord[0],mp_coord[1],mp_coord[2],current_time,0.0,0.0,0.0) * self.vector_direction

                        else:
                            self.value = self.vector_direction * self.modulus

                    mp.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION,[self.value],self.model_part.ProcessInfo)
                    mp.SetValuesOnIntegrationPoints(KratosMPM.MP_ACCELERATION,[self.value],self.model_part.ProcessInfo)
                    #node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, self.value[0])  # Set the x-component body force field
                    #node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, self.value[1])  


