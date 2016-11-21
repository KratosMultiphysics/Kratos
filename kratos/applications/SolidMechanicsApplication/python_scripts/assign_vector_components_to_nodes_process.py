import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given set by deffault to true.
import sys

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorComponentsToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class AssignVectorComponentsToNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help": "This process assigns a vector value to a vector variable component by component",
             "model_part_name": "MODEL_PART_NAME",
             "mesh_id": 0,
             "variable_name": "VARIABLE_NAME",           
             "value": [0.0, 0.0, 0.0],
             "constrained":true,
             "interval": [0.0, "End"],
             "local_axes" : {}
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
                
        ##check if variable type is a vector
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if( type(self.var) != KratosMultiphysics.Array1DVariable3 ):
            raise Exception("Variable type is incorrect. Must be a three-component vector.")

        self.model         = Model
        self.model_part    = Model[self.settings["model_part_name"].GetString()]
        self.variable_name = self.settings["variable_name"].GetString()

        ###check component assignation
        self.AssignValueProcesses = []
        
        self.constraints = []
        for i in range(0, self.settings["value"].size() ):
            if( self.settings["value"][i].IsNull() ):
                self.constraints.append(False)
            else:
                self.constraints.append(True)

        #print(" constraints ", self.constraints)
        
        self.BuildComponentsProcesses()
               
                      
    def ExecuteInitialize(self):
        for process in self.AssignValueProcesses:
            process.ExecuteInitialize()
 

    def ExecuteInitializeSolutionStep(self):
        for process in self.AssignValueProcesses:
            process.ExecuteInitializeSolutionStep()       

                
    def ExecuteFinalizeSolutionStep(self):
        for process in self.AssignValueProcesses:
            process.ExecuteFinalizeSolutionStep()
            
             
            
            
    #
    def BuildComponentsProcesses(self):

        counter = 0
        for imposed in self.constraints:
            
            if( imposed ):
                params = KratosMultiphysics.Parameters("{}")           
                params.AddValue("model_part_name", self.settings["model_part_name"])
                params.AddValue("mesh_id", self.settings["mesh_id"])
               
                if( counter == 0 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_X")
                if( counter == 1 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Y")
                if( counter == 2 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Z")

                params.AddValue("interval",self.settings["interval"])
                params.AddValue("constrained", self.settings["constrained"])
                
                if( self.settings["value"][counter].IsNumber() ):
                    params.AddEmptyValue("value").SetDouble(self.settings["value"][counter].GetDouble())
                else:
                    params.AddEmptyValue("value").SetString(self.settings["value"][counter].GetString())

                params.AddValue("local_axes", self.settings["local_axes"])                

                import assign_scalar_to_nodes_process as assign_scalar_process
                
                self.AssignValueProcesses.append(assign_scalar_process.AssignScalarToNodesProcess(self.model, params))
                
            counter +=1
  
            
        
        
        
