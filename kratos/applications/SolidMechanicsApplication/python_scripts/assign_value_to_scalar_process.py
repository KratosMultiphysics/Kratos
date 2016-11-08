import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import sys

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignValueToScalarProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class AssignValueToScalarProcess(KratosMultiphysics.Process):
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
             "value": 0.0
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
                
        ##check if variable type is a scalar or a vector component
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if(type(self.var) != KratosMultiphysics.Array1DComponentVariable and type(self.var) != KratosMultiphysics.DoubleVariable):
            raise Exception("Variable type is incorrect. Must be a scalar or a component")

        self.model_part    = Model[self.settings["model_part_name"].GetString()]
        self.variable_name = self.settings["variable_name"].GetString()

        self.interval  = []
        self.interval.append(self.settings["interval"][0].GetDouble());
        if( self.settings["interval"][1].IsString() ):
            if( self.settings["interval"][1].GetString() == "End" ):
                self.interval.append(sys.float_info.max)
        elif( self.settings["interval"][1].IsDouble() or  self.settings["interval"][1].IsInt() ):
            self.interval.append(self.settings["interval"][1].GetDouble());

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])
      
        self.function_string = self.settings["time_function"].GetString()
        if( self.function_string == "constant" ):
            self.function_expression = "1";
        elif( self.function_string == "incremental" ):
            self.function_expression = "t"
        else:
            self.function_expression = self.function_string;
        if (sys.version_info > (3, 0)):
            self.compiled_function = compile(self.function_expression, '', 'eval', optimize=2)
        else:
            self.compiled_function = compile(self.function_expression, '', 'eval')
            
        if( self.interval[0] == 0.0 and self.interval[1] == 0.0 ):
            self.interval_string = "initial"
        else:
            self.interval_string = ""

        self.value = self.settings["value"].GetDouble();        

        self.AssignValueProcess = KratosMultiphysics.Process()
        self.FixDofsProcess     = KratosMultiphysics.Process()
        self.FreeDofsProcess    = KratosMultiphysics.Process()
        self.Started            = False

    def function(self, t):
        return eval(self.compiled_function)

       
    def ExecuteInitialize(self):              
        
        params = KratosMultiphysics.Parameters("{}")           
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("mesh_id", self.settings["mesh_id"])
        params.AddValue("variable_name", self.settings["variable_name"])
        
        if( self.interval_string != "initial" ):
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcess.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcess.append(free_dof_process)
                
        params.AddValue("value", self.settings["value"])    
        scalar_process = KratosSolid.AssignValueToScalarVariableProcess(self.model_part, params)
        self.AssignValueProcess.append(scalar_process)

        if( self.interval_string == "initial" ):
            self.AssignValueProcesses.Execute()


    def ExecuteInitializeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001;

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ): 
            if( self.Started == False ):
                self.FixDofsProcess.Execute()
                self.Started = True
            else:
                interval_time = self.model_part.ProcessInfo[KratosMultiphysics.INTERVAL_END_TIME]
                previous_time = self.model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.TIME]
                if(previous_time == interval_time):
                    self.FixDofsProcess.Execute()

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ):
            if( self.function_string == "constant" ):
                self.AssignValueProcess.Execute()
            else:        
                params = KratosMultiphysics.Parameters("{}")           
                params.AddValue("model_part_name", self.settings["model_part_name"])
                params.AddValue("mesh_id", self.settings["mesh_id"])
                params.AddValue("variable_name", self.settings["variable_name"])           
                params.AddValue("value", self.settings["value"]) 
                scalar_value = self.value * self.function(current_time)
                params.AddEmptyValue("value").SetDouble(scalar_value)

                scalar_process = KratosSolid.AssignValueToScalarVariableProcess(self.model_part, params)
                scalar_process.Execute()
                
    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time   = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                
        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001;

        if( (current_time + delta_time) > (self.interval[1] + tolerance) ):
            self.FreeDofsProcesses.Execute()
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, current_time)
            
            
            

        
        
        
