import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import math

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignValueAndDirectionToVectorProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class AssignValueAndDirectionToVectorProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "model_part_name": "MODEL_PART_NAME",
             "mesh_id": 0,
             "variable_name": "VARIABLE_NAME",
             "interval": [0.0, "End"],
             "time_function": "constant",
             "modulus" : 0.0,
             "direction": [0.0, 0.0, 0.0]
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
                
        ##check if variable type is a vector
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if( type(self.var) != KratosMultiphysics.Array1DVariable3 ):
            raise Exception("Variable type is incorrect. Must be a three-component vector.")


        self.model_part    = Model[self.settings["model_part_name"].GetString()]
        self.variable_name = self.settings["variable_name"].GetString()

        self.interval  = []
        self.interval.append(self.settings["interval"][0].GetDouble());
        if( self.settings["interval"][1].IsString() ):
            if( self.settings["interval"][1].GetString() == "End" ):
                self.interval.append(sys.float_info.max)
        elif( self.settings["interval"][1].IsDouble() ):
            self.interval.append(self.settings["interval"][1].GetDouble());
        
        self.function_string   = self.settings["time_function"].GetString()
        if( self.function_string == "constant" ):
            self.function_expression = "1";
        elif( self.function_string == "incremental" ):
            self.function_expression = "t"
        else:
            self.function_expression = self.function_string;
        self.compiled_function = compile(self.function_expression, '', 'eval', optimize=2)
            
        self.interval_string = "custom"
        if( self.interval[0] == 0.0 and self.interval[1] == 0.0 ):
            self.interval_string = "initial"

        self.AssignValueProcesses = []
        self.FixDofsProcesses     = []
        self.FreeDofsProcesses    = []

        modulus = self.settings["modulus"].GetDouble();

        direction   = []
        scalar_prod = 0 
        for i in range(0, self.settings["direction"].size() ):
            direction.append( self.settings["direction"][i].GetDouble() )
            scalar_prod = scalar_prod + direction[i]*direction[i]

        norm = math.sqrt(scalar_prod)

        self.value = []
        if( norm != 0.0 ):
            for j in direction:
                self.value.append( j*modulus/norm )
        else:
            for j in direction:
                self.value.append(0.0)


    def function(self, t):
        return eval(self.compiled_function)

       
    def ExecuteInitialize(self):
        
        for i in range(0, len(self.value) ):
            params = KratosMultiphysics.Parameters("{}")           
            params.AddValue("model_part_name", self.settings["model_part_name"])
            params.AddValue("mesh_id", self.settings["mesh_id"])
            if( i == 0 ):
                params.AddEmptyValue("variable_name").SetString(self.variable_name+"_X")
            if( i == 1 ):
                params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Y")
            if( i == 2 ):
                params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Z")

            params.AddEmptyValue("value").SetDouble(self.value[i])

            component_process = KratosSolid.AssignValueToScalarVariableProcess(self.model_part, params)
            self.AssignValueProcesses.append(component_process)
            
        if( self.interval_string == "initial" ):
            for process in self.AssignValueProcesses:
                process.Execute()

    def ExecuteInitializeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001;

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ):
            if( self.function_string == "constant" ):
                for process in self.AssignValueProcesses:
                    process.Execute()
            else:        
                for i in range(0, len(self.value) ):
                    params = KratosMultiphysics.Parameters("{}")           
                    params.AddValue("model_part_name", self.settings["model_part_name"])
                    params.AddValue("mesh_id", self.settings["mesh_id"])
                    if( i == 0 ):
                        params.AddEmptyValue("variable_name").SetString(self.variable_name+"_X")
                    if( i == 1 ):
                        params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Y")
                    if( i == 2 ):
                        params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Z")
                                        
                    component_value = self.value[i] * self.function(current_time)
                    params.AddEmptyValue("value").SetDouble(component_value)
                    
                    component_process = KratosSolid.AssignValueToScalarVariableProcess(self.model_part, params)
                    component_process.Execute()


            
            

        
        
        
