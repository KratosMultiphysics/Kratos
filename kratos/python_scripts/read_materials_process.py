import KratosMultiphysics  
import sys
        
def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReadMaterialsProcess(Model, settings["Parameters"])


class ReadMaterialsProcess(KratosMultiphysics.Process):
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
        KratosMultiphysics.Process.__init__(self) 
        default_settings = KratosMultiphysics.Parameters("""
            {
            "materials_filename" : "please specify the file to be opened"
            }
            """
        )

        settings.ValidateAndAssignDefaults(default_settings)
        self.Model = Model

        with open(settings["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())
        
        for i in range(materials["properties"].size()):
            self._AssignPropertyBlock(materials["properties"][i])
        
        print("finished reading materials")
        
    def _GetVariable(self, my_string):
        """Return the python object of a Variable named by the string argument.

        Examples:
        variable = self._GetVariable("VELOCITY")
        variable = self._GetVariable("KratosMultiphysics.VELOCITY")
        variable = self._GetVariable("SUBSCALE_PRESSURE")
        variable = self._GetVariable("FluidDynamicsApplication.SUBSCALE_PRESSURE")
        variable = self._GetVariable("KratosMultiphysics.FluidDynamicsApplication.SUBSCALE_PRESSURE")
        """
        splitted = my_string.split(".")

        if len(splitted) == 0:
            raise Exception("Something wrong. Trying to split the string " + my_string)
        if len(splitted) > 3:
            raise Exception("Something wrong. String " + my_string + " has too many arguments")

        return KratosMultiphysics.KratosGlobals.GetVariable(splitted[-1]) # This also checks if the application has been imported

    def _GetConstitutiveLaw(self, my_string):
        """Return the python object of a Constitutive Law named by the string argument.

        Example:
        constitutive_law = self._GetConstitutiveLaw('KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw')
        model_part.GetProperties(prop_id).SetValue(CONSTITUTIVE_LAW, constitutive_law)
        """
        splitted = my_string.split(".")

        if(len(splitted) == 0):
            raise Exception("something wrong. Trying to split the string " + my_string)
        if(len(splitted) == 1):
            raise Exception("Please also provide the name of the application of constitutive law " + my_string)
        if len(splitted) > 3:
            raise Exception("Something wrong. String " + my_string + " has too many arguments")

        constitutive_law_name = splitted[-1]
        application_name = splitted[-2]

        if application_name == "KratosMultiphysics":
            return getattr(KratosMultiphysics, constitutive_law_name) 
        else:
            # Check that application was imported in the main script
            KratosMultiphysics.CheckRegisteredApplications(application_name)
            application = __import__("Kratos" + application_name)
            
            return getattr(application, constitutive_law_name)

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
           constitutive_law = self._GetConstitutiveLaw( mat["constitutive_law"]["name"].GetString())(mat["constitutive_law"]["Variables"])
        else:
           constitutive_law = self._GetConstitutiveLaw( mat["constitutive_law"]["name"].GetString())()
           
        prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law)
        
        # Add / override the values of material parameters in the properties
        for key, value in mat["Variables"].items():
            var = self._GetVariable(key)
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

            input_var = self._GetVariable(table["input_variable"].GetString())
            output_var = self._GetVariable(table["output_variable"].GetString())

            new_table = KratosMultiphysics.PiecewiseLinearTable()

            for i in range(table["data"].size()):
                new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())

            prop.SetTable(input_var,output_var,new_table)

'''
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
            data = materials["properties"][i]
    
            interpolation_required = False
            for key, value in data["Material"]["Variables"].items():
                if value.IsSubParameter():
                    if value.Has("@table"):
                        interpolation_required = True
                        break
                    else:
                        raise ValueError("Variable Dict is not valid!")

            if interpolation_required:
                self._AssignPropertyBlockInterpolated(data)
            else:
                self._AssignPropertyBlock(data)
        
        print("finished reading materials")
        
    def _GetItemFromModule(self, my_string):
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

        self._SetConstitutiveLaw(prop, mat)
        
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
            input_var = self._GetItemFromModule(table["input_variable"].GetString())
            output_var = self._GetItemFromModule(table["output_variable"].GetString())
            self._SetTable(prop, table, input_var, output_var)       
            table_name = key

    def _SetConstitutiveLaw(self, prop, mat):
        # Set the CONSTITUTIVE_LAW for the current properties.
        if "Variables" in mat["constitutive_law"].keys(): #pass the list of variables when constructing the constitutive law
           constitutive_law = self._GetItemFromModule( mat["constitutive_law"]["name"].GetString())(mat["constitutive_law"]["Variables"])
        else:
           constitutive_law = self._GetItemFromModule( mat["constitutive_law"]["name"].GetString())()
           
        prop.SetValue(CONSTITUTIVE_LAW, constitutive_law)

    def _SetTable(self, prop, table, input_var, output_var):
        new_table = PiecewiseLinearTable()

        for i in range(table["data"].size()):
            new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())

        prop.SetTable(input_var,output_var,new_table)
             
    def _AssignPropertyBlockInterpolated(self, data):        
        """Set constitutive law and material properties and assign to elements and conditions.
        
        Some variables are interpolated from Tables

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
                    "YOUNG_MODULUS" : {"@table : "E_Values"},
                    "POISSON_RATIO" : 0.3,
                    "RESIDUAL_VECTOR" : [1.5,0.3,-2.58],
                    "LOCAL_INERTIA_TENSOR" : [[0.27,0.0],[0.0,0.27]]
                },
                "Tables" : {
                    ""E_Values" : {
                        "input_variable" : "TEMPERATURE",
                        "output_variable" : "YOUNG_MODULUS",
                        "input_variable_location" : "geom_entity",
                        "data" : [
                            [0.0,  100.0],
                            [20.0, 90.0],
                            [30.0, 85.0],
                            [35.0, 80.0]
                        ]
                    }
                }
            }
        }
        """
        # Get the properties for the specified model part.
        model_part = self.Model[data["model_part_name"].GetString()]
        mesh_id = 0
        mat = data["Material"]

        # Assign the properties to the model part's elements and conditions.
        for elem in model_part.Elements:
            current_number_props = model_part.NumberOfProperties()
            elem_props = model_part.GetProperties(current_number_props + 1, mesh_id)

            elem.Properties = elem_props

            self._AssignInterpolatedProps(mat, elem, elem_props)

        for cond in model_part.Conditions:
            current_number_props = model_part.NumberOfProperties()
            cond_props = model_part.GetProperties(current_number_props + 1, mesh_id)

            cond.Properties = cond_props

            self._AssignInterpolatedProps(mat, cond, cond_props)

    def _AssignInterpolatedProps(self, mat, geom_entity, prop):
        self._SetConstitutiveLaw(prop, mat)
        
        # Add / override tables in the properties
        tables_for_interpolation = {}
        for key, table in mat["Tables"].items():
            input_var = self._GetItemFromModule(table["input_variable"].GetString())
            output_var = self._GetItemFromModule(table["output_variable"].GetString())
            self._SetTable(prop, table, input_var, output_var)  
            table_name = key

            # save information abt the tables used for interpolation
            if table.Has("input_variable_location"):
                input_var_location = table["input_variable_location"].GetString()
                tables_for_interpolation[table_name] = [input_var,output_var,input_var_location]

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
                prop.SetValue( var, string_val )
            elif value.IsMatrix():
                prop.SetValue( var, value.GetMatrix() )
            elif value.IsVector():
                prop.SetValue( var, value.GetVector() )
            elif value.IsSubParameter():
                if value.Has("@table"):
                    table_name = value["@table"].GetString()
                    if table_name not in tables_for_interpolation.keys():
                        raise KeyError("Table \"" + table_name + "\" needed for interpolation is not defined!")
                    table_info = tables_for_interpolation[table_name]
                    interpolated_value = self._ComputeInterpolatedValue(geom_entity, prop, table_info, var)
                    prop.SetValue( var, interpolated_value )
                else:
                    raise ValueError("Variable keyword is not valid!")
            else:
                raise ValueError("Type of value is not available")

    def _ComputeInterpolatedValue(self, geom_entity, prop, table_info, dict_key):
        input_var = table_info[0]
        input_var_location = table_info[2]
        output_var = table_info[1]

        if dict_key != output_var:
            raise Exception("Variable mismatch!")

        table_to_interpolate = prop.GetTable(input_var, output_var)

        if input_var_location == "geom_entity":
            if geom_entity.Has(input_var): # Values in Geom Entites are saved as Non-historical values (model_part_io.cpp)
                input_value = geom_entity.GetValue(input_var)
            else:
                err_msg  = "Geometric Entity # " + str(geom_entity.Id)
                err_msg += " does not have " + str(input_var)
                raise ValueError(err_msg)
        elif input_var_location == "nodes":
            nodes = geom_entity.GetNodes()
            input_value = 0
            for node in nodes:
                if node.SolutionStepsDataHas(input_var): # Values in Nodes are saved as Historical values (model_part_io.cpp)
                    input_value += node.GetSolutionStepValue(input_var)
                else:
                    err_msg  = "Node # " + str(node.Id)
                    err_msg += " does not have " + str(input_var)
                    raise ValueError(err_msg)
                input_value /= len(nodes)
        else:
            raise Exception("Type of input_var_location \"" + input_var_location + "\" is not valid!")

        interpolated_value = table_to_interpolate.GetValue(input_value) # interpolate the value from the table

        return interpolated_value
'''