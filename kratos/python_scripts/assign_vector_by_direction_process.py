# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import assign_scalar_variable_process

import math

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorByDirectionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorByDirectionProcess(KratosMultiphysics.Process):
    """This process sets a variable a certain scalar value in a given direction, for all the nodes belonging to a submodelpart. Uses assign_scalar_variable_to_conditions_process for each component

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                 : "This process sets a variable a certain scalar value in a given direction, for all the nodes belonging to a submodelpart. Uses assign_scalar_variable_to_conditions_process for each component",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "modulus"              : 0.0,
            "constrained"          : true,
            "direction"            : [1.0, 0.0, 0.0],
            "local_axes"           : {}
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if settings.Has("modulus"):
            if settings["modulus"].IsString():
                default_settings["modulus"].SetString("0.0")
        else:
            raise RuntimeError("Please specify the modulus of the vector")

        if settings.Has("direction"):
            if settings["direction"].IsString():
                default_settings["direction"].SetString("Automatic")
        else:
            raise RuntimeError("Please specify the direction of the vector")

        # Detect "End" as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" + settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # Construct the component by component parameter objects
        x_params = KratosMultiphysics.Parameters("{}")
        y_params = KratosMultiphysics.Parameters("{}")
        z_params = KratosMultiphysics.Parameters("{}")

        list_params = [x_params, y_params, z_params]
        for i_dir, var_string in enumerate(["_X", "_Y", "_Z"]):
            list_params[i_dir].AddValue("model_part_name",settings["model_part_name"])
            list_params[i_dir].AddValue("mesh_id",settings["mesh_id"])
            list_params[i_dir].AddValue("constrained",settings["constrained"])
            list_params[i_dir].AddValue("interval",settings["interval"])
            list_params[i_dir].AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + var_string)
            list_params[i_dir].AddValue("local_axes",settings["local_axes"])

        # "Automatic" direction: get the inwards direction
        all_numeric = True
        if settings["direction"].IsString():
            if settings["direction"].GetString() == "automatic_inwards_normal" or settings["direction"].GetString() == "automatic_outwards_normal":
                # Compute the condition normals
                KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

                # Compute the average conditions normal in the submodelpart of interest
                avg_normal = KratosMultiphysics.VariableUtils().SumConditionVectorVariable(KratosMultiphysics.NORMAL, self.model_part)
                avg_normal_norm = math.sqrt(pow(avg_normal[0],2) + pow(avg_normal[1],2) + pow(avg_normal[2],2))
                if avg_normal_norm < 1.0e-12:
                    raise Exception("Direction norm is close to 0 in AssignVectorByDirectionProcess.")

                unit_direction = KratosMultiphysics.Vector(3)
                unit_direction = (1.0/avg_normal_norm) * avg_normal

                # Note that the NormalCalculationUtils().CalculateOnSimplex gives the outwards normal vector
                if settings["direction"].GetString() == "automatic_inwards_normal":
                    unit_direction = (-1)*unit_direction
        # Direction is given as a vector
        elif settings["direction"].IsArray():
            unit_direction = [0.0,0.0,0.0]
            direction_norm = 0.0
            for i in range(0,3):
                if settings["direction"][i].IsNumber():
                    unit_direction[i] = settings["direction"][i].GetDouble()
                    direction_norm += pow(unit_direction[i],2)
                else:
                    function_string = settings["direction"][i].GetString()
                    unit_direction[i] = function_string
                    all_numeric = False

            # Normalize direction
            if all_numeric:
                direction_norm = math.sqrt(direction_norm)
                if direction_norm < 1.0e-12:
                    raise Exception("Direction norm is close to 0 in AssignVectorByDirectionProcess.")
                for i in range(0,3):
                    unit_direction[i] = unit_direction[i]/direction_norm

        # Set the remainding parameters
        if settings["modulus"].IsNumber():
            modulus = settings["modulus"].GetDouble()
            if all_numeric:
                for i_dir in range(3):
                    list_params[i_dir].AddEmptyValue("value").SetDouble(modulus * unit_direction[i_dir])
            else:
                for i_dir in range(3):
                    list_params[i_dir].AddEmptyValue("value").SetString("("+str(unit_direction[i_dir])+")*("+str(modulus)+")")
        elif settings["modulus"].IsString():
            # The concatenated string is: "direction[i])*(f(x,y,z,t)"
            modulus = settings["modulus"].GetString()
            for i_dir in range(3):
                list_params[i_dir].AddEmptyValue("value").SetString("("+str(unit_direction[i_dir])+")*("+modulus+")")


        # Construct a AssignScalarToNodesProcess for each component
        self.aux_processes = []
        for i_dir in range(3):
            self.aux_processes.append( assign_scalar_variable_process.AssignScalarVariableProcess(Model, list_params[i_dir]) )

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed in before initialize the solution loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()
