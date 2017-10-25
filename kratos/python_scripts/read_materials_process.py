from KratosMultiphysics import * 
import importlib
        
def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReadMaterialsProcess(Model, settings["Parameters"])


class ReadMaterialsProcess(Process):
    def __init__(self, Model, settings):
        """Read constitutive law and material properties from a json file and assign them to elements and conditions.

        Arguments:
        Model -- a dictionary of model parts to which properties may be assigned.
        settings -- Kratos parameters object specifying the name of the json file containing property information.

        Example:
        params = Parameters('{"materials_filename" : "materials.json"}')
        read_materials_process = ReadMaterialsProcess(params)

        The json file must include a list of properties:
        {
            "properties" : []
        }

        See _AssignPropertyBlock for detail on how properties are imported.
        """
        Process.__init__(self) 
        default_settings = Parameters("""
            {
            "materials_filename" : "please specify the file to be opened"
            }
            """
            )
            
        settings.ValidateAndAssignDefaults(default_settings)
        self.Model = Model

        parameter_file = open(settings["materials_filename"].GetString(), 'r')
        materials = Parameters(parameter_file.read())
        
        for i in range(materials["properties"].size()):
            self._AssignPropertyBlock(materials["properties"][i])
        
        print("finished reading materials")
        
    def _GetItemFromModule(self,my_string):
        """Return the python object named by the string argument.

        Example:
        constitutive_law = self._GetItemFromModule('KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw')
        model_part.GetProperties(prop_id).SetValue(CONSTITUTIVE_LAW, constitutive_law)
        """
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
        """Set constitutive law and material properties and assign to elements and conditions.

        Arguments:
        data -- a dictionary or json object defining properties for a model part.

        Example:
        data = {
            "model_part_name" : "Plate",
            "properties_id" : 1,
            "Material" : {
                "constitutive_law" : {
                    "name" : "KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw"
                },
                "Variables" : {
                    "YOUNG_MODULUS" : 200e9,
                    "POISSON_RATIO" : 0.3,
                    "RESIDUAL_VECTOR" : [1.5,0.3,-2.58],
                    "LOCAL_INERTIA_TENSOR" : [[0.27,0.0],[0.0,0.27]]
                },
                "Tables" : {}
            }
        }
        """
        # Get the properties for the specified model part.
        model_part = self.Model[data["model_part_name"].GetString()]
        property_id = data["properties_id"].GetInt()
        mesh_id = 0
        prop = model_part.GetProperties(property_id, mesh_id)
        
        # Assign the properties to the model part's elements and conditions.
        for elem in model_part.Elements:
            elem.Properties = prop
            
        for cond in model_part.Conditions:
            cond.Properties = prop
        
        mat = data["Material"]

        # Set the CONSTITUTIVE_LAW for the current properties.
        if "Variables" in mat["constitutive_law"].keys(): #pass the list of variables when constructing the constitutive law
           constitutive_law = self._GetItemFromModule( mat["constitutive_law"]["name"].GetString())(mat["constitutive_law"]["Variables"])
        else:
           constitutive_law = self._GetItemFromModule( mat["constitutive_law"]["name"].GetString())()
           
        prop.SetValue(CONSTITUTIVE_LAW, constitutive_law)
        
        # Add / override the values of material parameters in the properties
        for key, value in mat["Variables"].items():
            var = self._GetItemFromModule(key)
            if value.IsDouble():
                prop.SetValue( var, value.GetDouble() )
            elif value.IsInt():
                prop.SetValue( var, value.GetInt() )
            elif value.IsBool():
                prop.SetValue( var, value.GetBool() )
            elif value.IsString():
                prop.SetValue( var, value.GetString() )
            elif value.IsMatrix():
                prop.SetValue( var, value.GetMatrix() )
            elif value.IsVector():
                prop.SetValue( var, value.GetVector() )
            else:
                raise ValueError("Type of value is not available")

        # Add / override tables in the properties
        for key, table in mat["Tables"].items():
            table_name = key

            input_var = self._GetItemFromModule(table["input_variable"].GetString())
            output_var = self._GetItemFromModule(table["output_variable"].GetString())

            new_table = PiecewiseLinearTable()

            for i in range(table["data"].size()):
                new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())

            prop.SetTable(input_var,output_var,new_table)
