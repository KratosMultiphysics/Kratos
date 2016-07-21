from KratosMultiphysics import * 
from KratosMultiphysics.SolidMechanicsApplication import *
        
def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReadMaterialsProcess(Model, settings["Parameters"])

    

##all the processes python processes should be derived from "python_process"
class ReadMaterialsProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self) 
        default_settings = Parameters("""
            {
            "materials_filename" : "please specify the file to be opened"
            }
            """
            )
            
        settings.ValidateAndAssignDefaults(default_settings)
        self.Model = Model
        
        #TODO: change to use the KratosParameters once dictionary iterators are exported
        import json
        with open(settings["materials_filename"].GetString()) as data_file:    
            materials = json.load(data_file)

        for prop in materials["properties"]:
            self._AssignPropertyBlock(prop)
        
        print("finished reading materials")
            
            
    def _AssignPropertyBlock(self, data):
        print(data["model_part_name"])
        model_part = self.Model[data["model_part_name"]]
        prop = model_part.Properties[ data["properties_id"] ]
                
        ####################################
        #assign property to the list of elements and conditions in the model part
        for elem in model_part.Elements:
            elem.Properties = prop
            
        for cond in model_part.Conditions:
            cond.Properties = prop
        
        
        ####################################
        ## read material data
        
        mat = data["Material"]

        #read constitutive law and assign it to prop
        if "Variables" in mat["constitutive_law"].keys(): #pass the list of variables when constructing the constitutive law
           constitutive_law = eval( mat["constitutive_law"]["name"])(mat["constitutive_law"]["Variables"])
        else:
           constitutive_law = eval( mat["constitutive_law"]["name"])()
           
        prop.SetValue(CONSTITUTIVE_LAW, constitutive_law)
        
        #read variables 
        for key, value in mat["Variables"].items():
            prop.SetValue(eval(key), value)
            #print(prop)

        #read table
        #for table in mat["Tables"]:
            #new_table = PiecewiseLinearTable(eval(table["input_variable"]), eval(table["output_variable"]))
            #for i in range(table["data"].size()):
                #new_table.AddRow(new_table["data"][i][0], new_table["data"][i][1])
            #prop.AddTable(new_table)
      
        
                
        
