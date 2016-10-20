import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.
import sys

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignValueToVectorComponentsProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class AssignValueToVectorComponentsProcess(KratosMultiphysics.Process):
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
             "imposed_components": [false, false, false],
             "value": [0.0, 0.0, 0.0]
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
           
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])

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

        self.constraints = self.settings["imposed_components"]

        self.AssignValueProcesses = []
        self.FixDofsProcesses     = []
        self.FreeDofsProcesses    = []
        self.Started              = False
 

        for i in range(0, self.constraints.size() ):

            if( self.settings["imposed_components"][i].GetBool() ):
                params = KratosMultiphysics.Parameters("{}")           
                params.AddValue("model_part_name", self.settings["model_part_name"])
                params.AddValue("mesh_id", self.settings["mesh_id"])
                if( i == 0 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_X")
                if( i == 1 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Y")
                if( i == 2 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Z")

                if( self.interval_string != "initial" ):
                    fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
                    self.FixDofsProcesses.append(fix_dof_process)
                    free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
                    self.FreeDofsProcesses.append(free_dof_process)
                
                params.AddEmptyValue("value").SetDouble(self.settings["value"][i].GetDouble())
                component_process = KratosSolid.AssignValueToScalarVariableProcess(self.model_part, params)
                self.AssignValueProcesses.append(component_process)

        

        # in dynamic problems derivated variables can be fixed
        self.fix_derivated_variable = False
        if( self.variable_name == "ACCELERATION" or self.variable_name == "VELOCITY" ):
            self.derivated_variable_name = "DISPLACEMENT"
            self.fix_derivated_variable = True

        if( self.variable_name == "ANGULAR_ACCELERATION" or self.variable_name == "ANGULAR_VELOCITY" ):
            self.derivated_variable_name = "ROTATION"
            self.fix_derivated_variable = True
            
        if( self.fix_derivated_variable ):

            for i in range(0, self.constraints.size() ):

                if( self.settings["imposed_components"][i].GetBool() ):
                    params = KratosMultiphysics.Parameters("{}")           
                    params.AddValue("model_part_name", self.settings["model_part_name"])
                    params.AddValue("mesh_id", self.settings["mesh_id"])
                    if( i == 0 ):
                        params.AddEmptyValue("variable_name").SetString(self.derivated_variable_name+"_X")
                    if( i == 1 ):
                        params.AddEmptyValue("variable_name").SetString(self.derivated_variable_name+"_Y")
                    if( i == 2 ):
                        params.AddEmptyValue("variable_name").SetString(self.derivated_variable_name+"_Z")

                    if( self.interval_string != "initial" ):
                        fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
                        self.FixDofsProcesses.append(fix_dof_process)
                        free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
                        self.FreeDofsProcesses.append(free_dof_process)
                    

    def function(self, t):
        return eval(self.compiled_function)

       
    def ExecuteInitialize(self):              
        
        if( self.interval_string == "initial" ):
            for process in self.AssignValueProcesses:
                process.Execute()

    def ExecuteInitializeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001;

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ):
            if( self.Started == False ):
                for process in self.FixDofsProcesses:
                    process.Execute()
                self.Started = True
            else:
                interval_time = self.model_part.ProcessInfo[KratosMultiphysics.INTERVAL_END_TIME]
                previous_time = self.model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.TIME]
                if(previous_time == interval_time):
                    for process in self.FixDofsProcesses:
                        process.Execute()

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ):
            if( self.function_string == "constant" ):
                for process in self.AssignValueProcesses:
                    process.Execute()
            else:        
                for i in range(0, self.constraints.size() ):
                    if( self.settings["imposed_components"][i].GetBool() ):
                        params = KratosMultiphysics.Parameters("{}")           
                        params.AddValue("model_part_name",self.settings["model_part_name"])
                        params.AddValue("mesh_id",self.settings["mesh_id"])
                        if( i == 0 ):
                            params.AddEmptyValue("variable_name").SetString(self.variable_name+"_X")
                        if( i == 1 ):
                            params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Y")
                        if( i == 2 ):
                            params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Z")
                                
                        component_value = self.settings["value"][i].GetDouble() * self.function(current_time)

                        #print("::ASSIGN::",self.model_part.Name," ",params["variable_name"].GetString(),": ",component_value)
                        params.AddEmptyValue("value").SetDouble(component_value)
                        component_process = KratosSolid.AssignValueToScalarVariableProcess(self.model_part, params)
                        component_process.Execute()


                
    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001;

        if( (current_time + delta_time) > (self.interval[1] + tolerance) ):
            for process in self.FreeDofsProcesses:
                process.Execute()
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, current_time)
            
            
            
            

        
        
        
