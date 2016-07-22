from KratosMultiphysics import * 
import importlib
#from KratosMultiphysics.SolidMechanicsApplication import *
        
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
        
    def _GetItemFromModule(self,my_string):
        splitted = my_string.split(".")
        if(len(splitted) == 0):
            raise Exception("something wrong. Trying to split the string "+my_string)
        if(len(splitted) == 1):
            return eval(my_string)
        else:
            module_name = ""
            for i in range(len(splitted)-1):
                module_name += splitted[i] 
                if i != len(splitted)-2:
                    module_name += "."

            module = importlib.import_module(module_name)
            return getattr(module,splitted[-1]) 
            
            
            
    def _AssignPropertyBlock(self, data):
        model_part = self.Model[data["model_part_name"]]
        property_id = data["properties_id"]
        mesh_id = 0
        prop = model_part.GetProperties(property_id, mesh_id)
        
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
           constitutive_law = self._GetItemFromModule( mat["constitutive_law"]["name"])(mat["constitutive_law"]["Variables"])
        else:
           constitutive_law = self._GetItemFromModule( mat["constitutive_law"]["name"])()
           
        prop.SetValue(CONSTITUTIVE_LAW, constitutive_law)
        
        #read variables 
        for key, value in mat["Variables"].items():
            var = self._GetItemFromModule(key)
            prop.SetValue( var, value)

        #read table
        for key, table in mat["Tables"].items():
            table_name = key

            input_var = self._GetItemFromModule(table["input_variable"])
            output_var = self._GetItemFromModule(table["output_variable"])

            new_table = PiecewiseLinearTable()
            for i in range(len(table["data"])):
                new_table.AddRow(table["data"][i][0], table["data"][i][1])
            prop.SetTable(input_var,output_var,new_table)
      
        
                
        
