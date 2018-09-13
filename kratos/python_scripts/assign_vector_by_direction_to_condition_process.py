import math
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorByDirectionToConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorByDirectionToConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                 : "This process sets a variable a certain scalar value in a given direction, for all the conditions belonging to a submodelpart. Uses assign_scalar_variable_to_conditions_process for each component",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
                "interval"             : [0.0, 1e30],
                "modulus"              : 1.0,
                "direction"            : [1.0, 0.0, 0.0],
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

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # Construct the component by component parameter objects
        x_params = KratosMultiphysics.Parameters("{}")
        y_params = KratosMultiphysics.Parameters("{}")
        z_params = KratosMultiphysics.Parameters("{}")

        # Component X
        x_params.AddValue("model_part_name",settings["model_part_name"])
        x_params.AddValue("mesh_id",settings["mesh_id"])
        x_params.AddValue("interval",settings["interval"])
        x_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_X")
        x_params.AddValue("local_axes",settings["local_axes"])

        # Component Y
        y_params.AddValue("model_part_name",settings["model_part_name"])
        y_params.AddValue("mesh_id",settings["mesh_id"])
        y_params.AddValue("interval",settings["interval"])
        y_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Y")
        y_params.AddValue("local_axes",settings["local_axes"])

        # Component Z
        z_params.AddValue("model_part_name",settings["model_part_name"])
        z_params.AddValue("mesh_id",settings["mesh_id"])
        z_params.AddValue("interval",settings["interval"])
        z_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Z")
        z_params.AddValue("local_axes",settings["local_axes"])

        # "Automatic" direction: get the inwards direction
        if(settings["direction"].IsString()):
            if ((settings["direction"].GetString() == "automatic_inwards_normal") or (settings["direction"].GetString() == "automatic_outwards_normal")):
                # Compute the condition normals
                KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

                # Compute the average conditions normal in the submodelpart of interest
                avg_normal = KratosMultiphysics.VariableUtils().SumConditionVectorVariable(KratosMultiphysics.NORMAL, self.model_part)
                avg_normal_norm = math.sqrt(pow(avg_normal[0],2) + pow(avg_normal[1],2) + pow(avg_normal[2],2))
                if(avg_normal_norm < 1e-6):
                    raise Exception("Direction norm is close to 0 in AssignVectorByDirectionToConditionProcess.")

                unit_direction = KratosMultiphysics.Vector(3)
                unit_direction = (1/avg_normal_norm)*avg_normal

                # Note that the NormalCalculationUtils().CalculateOnSimplex gives the outwards normal vector
                if (settings["direction"].GetString() == "automatic_inwards_normal"):
                    unit_direction = (-1)*unit_direction

        # Direction is given as a vector
        elif(settings["direction"].IsArray()):
            # Normalize direction
            direction_norm = math.sqrt(pow(settings["direction"][0].GetDouble(),2) +
                                       pow(settings["direction"][1].GetDouble(),2) +
                                       pow(settings["direction"][2].GetDouble(),2))
            if(direction_norm < 1e-6):
                raise Exception("Direction norm is close to 0 in AssignVectorByDirectionToConditionProcess.")

            unit_direction = []
            for i in range(0,3):
                unit_direction.append(settings["direction"][i].GetDouble()/direction_norm)


        # Set the remainding parameters
        if(settings["modulus"].IsNumber()):
            modulus = settings["modulus"].GetDouble()
            x_params.AddEmptyValue("value").SetDouble(modulus*unit_direction[0])
            y_params.AddEmptyValue("value").SetDouble(modulus*unit_direction[1])
            z_params.AddEmptyValue("value").SetDouble(modulus*unit_direction[2])
        elif(settings["modulus"].IsString()):
            # The concatenated string is: "direction[i])*(f(x,y,z,t)"
            modulus = settings["modulus"].GetString()
            x_params.AddEmptyValue("value").SetString("("+str(unit_direction[0])+")*("+modulus+")")
            y_params.AddEmptyValue("value").SetString("("+str(unit_direction[1])+")*("+modulus+")")
            z_params.AddEmptyValue("value").SetString("("+str(unit_direction[2])+")*("+modulus+")")

        # Construct a AssignScalarToNodesProcess for each component
        import assign_scalar_variable_to_conditions_process

        self.aux_processes = []
        self.aux_processes.append( assign_scalar_variable_to_conditions_process.AssignScalarVariableToConditionsProcess(Model, x_params) )
        self.aux_processes.append( assign_scalar_variable_to_conditions_process.AssignScalarVariableToConditionsProcess(Model, y_params) )
        self.aux_processes.append( assign_scalar_variable_to_conditions_process.AssignScalarVariableToConditionsProcess(Model, z_params) )

    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()
