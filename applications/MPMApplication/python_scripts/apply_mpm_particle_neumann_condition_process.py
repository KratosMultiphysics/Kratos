import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMParticleNeumannConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMParticleNeumannConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"                 : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "material_points_per_condition"   : 0,
                "variable_name"                   : "PLEASE_SPECIFY_LOADING_CONDITION",
                "modulus"                         : 1.0,
                "constrained"                     : "fixed",
                "direction"                       : [0.0, 0.0, 0.0],
                "interval"                        : [0.0, 1e30],
                "option"                          : "",
                "local_axes"                      : {}
            }  """ )

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if settings.Has("modulus"):
            if settings["modulus"].IsString():
                default_parameters["modulus"].SetString("0.0")

        if settings.Has("direction"):
            if settings["direction"].IsString():
                default_parameters["direction"].SetString("Automatic")

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        context_string = type(self).__name__
        old_name = 'particles_per_condition'
        new_name = 'material_points_per_condition'
        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.model_part_name = settings["model_part_name"].GetString()
        self.model = Model
        self.material_points_per_condition = settings["material_points_per_condition"].GetInt()
        self.is_neumann_boundary = True
        self.option = settings["option"].GetString()

        # check constraint
        self.constrained = settings["constrained"].GetString()
        if (self.constrained == "fixed"):
            self.normal_following_load = False
        elif (self.constrained == "normal"):
            self.normal_following_load = True
        else:
            err_msg =  "The requested type of constrain: \"" + self.constrained + "\" is not available!\n"
            err_msg += "Available options are: \"fixed\" and \"normal\"."
            raise Exception(err_msg)

        # get variable imposed and check
        variable_name = settings["variable_name"].GetString()
        variable_name_list = ["POINT_LOAD","LINE_LOAD","SURFACE_LOAD"]
        if(variable_name in variable_name_list):
            self.variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        else:
            err_msg =  "The given variable \"" + variable_name + "\" is not available to be imposed with this process.\n"
            err_msg += "Available options are: " + ", ".join(variable_name_list)
            raise Exception(err_msg)


        self.vector_direction = KratosMultiphysics.Vector(3)
        self.aux_function_direction = ["0.0","0.0","0.0"]
        self.value_direction_is_numeric = [False, False, False]
        for i in range(3):
            if settings["direction"][i].IsNumber():
                self.value_direction_is_numeric[i] = True
                self.vector_direction[i] = settings["direction"][i].GetDouble()
            else:
                self.function_string_direction = settings["direction"][i].GetString()
                self.aux_function_direction[i] = KratosMultiphysics.GenericFunctionUtility(self.function_string_direction, settings["local_axes"])


        self.value_is_numeric = False
        self.value = KratosMultiphysics.Vector(3)
        self.aux_function = "0.0"

        if settings["modulus"].IsNumber():
            self.value_is_numeric = True
            self.modulus = settings["modulus"].GetDouble()
            self.value = self.vector_direction * self.modulus
        else:
            self.function_string = settings["modulus"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])


        self.modified_normal = False

        # Set Flag BOUNDARY and variables MATERIAL_POINTS_PER_CONDITION
        if self.material_points_per_condition >= 0:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, self.model_part.Nodes)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.BOUNDARY, True)
                condition.Set(KratosMultiphysics.MARKER, self.normal_following_load)
                condition.Set(KratosMultiphysics.MODIFIED, self.modified_normal)
                condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, self.material_points_per_condition)
                condition.SetValue(KratosMPM.MPC_IS_NEUMANN, self.is_neumann_boundary)
                condition.SetValue(self.variable, self.value)
        else:
            err_msg = '\n::[ApplyMPMParticleNeumannConditionProcess]:: W-A-R-N-I-N-G: You have specified invalid "material_points_per_condition", '
            err_msg += 'or assigned negative values. \nPlease assign: "material_points_per_condition" > 0 or = 0 (for automatic value)!\n'
            raise Exception(err_msg)

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed in before initialize the solution step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        for mpc in self.model_part.Conditions:
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

            mpc_coord = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]


            if self.interval.IsInInterval(current_time):

                self.step_is_active = True

                for i in range(3):
                    if not self.value_direction_is_numeric[i]:
                        if self.aux_function_direction[i].DependsOnSpace() == False: #depends on time only
                            self.vector_direction[i] = self.aux_function_direction[i].CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                        else: #most general case - space varying function (possibly also time varying)
                            self.vector_direction[i]= self.aux_function_direction[i].CallFunction(mpc_coord[0],mpc_coord[1],mpc_coord[2],current_time,0.0,0.0,0.0)


                if not self.value_is_numeric:
                    if self.aux_function.DependsOnSpace() == False: #depends on time only
                        self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0) * self.vector_direction
                    else: #most general case - space varying function (possibly also time varying)
                        self.value = self.aux_function.CallFunction(mpc_coord[0],mpc_coord[1],mpc_coord[2],current_time,0.0,0.0,0.0) * self.vector_direction
                else:
                    self.value = self.vector_direction * self.modulus


                mpc.SetValuesOnIntegrationPoints(KratosMPM.POINT_LOAD,[self.value],self.model_part.ProcessInfo)

