#from cmath import pi, sin
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
#from KratosMultiphysics.ParticleMechanicsApplication import assign_scalar_initial_condition_to_particle_process
#import math

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInitialConditionToParticleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignInitialConditionToParticleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
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
        #self.variable = KratosParticle.MP_VELOCITY
        if not isinstance(self.variable, KratosMultiphysics.Array1DVariable3) and not isinstance(self.variable, KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        settings.ValidateAndAssignDefaults(default_settings)

        # Get updated model_part
        self.model = Model
        model_part_name = settings["model_part_name"].GetString()
        if (model_part_name.startswith('Initial_MPM_Material.')):
            model_part_name = model_part_name.replace('Initial_MPM_Material.','')
            self.mpm_material_model_part_name = "MPM_Material." + model_part_name
        else:
            self.mpm_material_model_part_name = model_part_name

        # The actual initial velocity application occurs after the submodelpart is
        # transferred from the initial MPM material to the MPM material in the particle
        # generator utility. Therefore we change the prefix from initial MPM material
        # to MPM material.


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
        self.velocity = KratosMultiphysics.Vector(3)
        self.aux_function = "0.0"

        if settings["modulus"].IsNumber():
            self.value_is_numeric = True
            self.modulus = settings["modulus"].GetDouble()
            self.value = self.vector_direction * self.modulus
        else:
            self.function_string = settings["modulus"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])


    def ExecuteBeforeSolutionLoop(self):
        # Assign velocity to MP after solver.Initialize() - only apply once at the beginning!
        model_part = self.model[self.mpm_material_model_part_name]
        current_time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        # the model part is identified here, AFTER it has been transferred to the MPM_material part!
        aux=0;
        for mp in model_part.Elements:

            print(mp)

            mp_coord = mp.CalculateOnIntegrationPoints(KratosParticle.MP_COORD,model_part.ProcessInfo)[0]

            for i in range(3):
                if not self.value_direction_is_numeric[i]:
                    if self.aux_function_direction[i].DependsOnSpace() == False: #depends on time only
                        self.vector_direction[i] = self.aux_function_direction[i].CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                    else: #most general case - space varying function (possibly also time varying)
                        self.vector_direction[i]= self.aux_function_direction[i].CallFunction(mp_coord[0],mp_coord[1],mp_coord[2],current_time,0.0,0.0,0.0)


            if not self.value_is_numeric:
                if self.aux_function.DependsOnSpace() == False: #depends on time only
                    self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0) * self.vector_direction
                else: #most general case - space varying function (possibly also time varying)
                    self.value = self.aux_function.CallFunction(mp_coord[0],mp_coord[1],mp_coord[2],current_time,0.0,0.0,0.0) * self.vector_direction

            else:
                self.value = self.vector_direction * self.modulus

            mp.SetValuesOnIntegrationPoints(self.variable,[self.value],model_part.ProcessInfo)
