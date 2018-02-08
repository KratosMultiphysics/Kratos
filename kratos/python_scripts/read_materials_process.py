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
            data = materials["properties"][i]

            if self.__interpolation_is_required(data):
                self._AssignPropertyBlockInterpolated(data)
            else:
                self._AssignPropertyBlock(data)

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
                    "RESIDUAL_VECTOR" : [1.5,{ "@table": "Table3" },-2.58],
                    "LOCAL_INERTIA_TENSOR" : [[0.27,{ "@table": "Table3" }],[0.0,0.27]]
                },
                "Tables" : {
                    "E_Values" : {
                        "input_variable" : "TEMPERATURE",
                        "output_variable" : "YOUNG_MODULUS",
                        "input_variable_location" : "geom_entity",
                        "data" : [
                            [0.0,  100.0],
                            [20.0, 90.0],
                            [30.0, 85.0],
                            [35.0, 80.0]
                        ]
                    },
                    "Table3": {
                        "input_variable": "AUX_INDEX",
                        "output_variable": "CAUCHY_STRESS_VECTOR",
                        "data": [
                            [ 0.0,   80.0 ],
                            [ 200.0, 190.0 ]
                        ],
                        "input_variable_location": "nodes"
                    }
                }
            }
        }
        """
        # Get the properties for the specified model part.
        model_part = self.Model[data["model_part_name"].GetString()]
        root_model_part = model_part.GetRootModelPart()
        mesh_id = 0
        mat = data["Material"]

        # Get the Constitutive Law
        constitutive_law = self._GetConstitutiveLaw( mat["constitutive_law"]["name"].GetString() )()

        # Set the tables to the ModelPart
        table_dict = {}
        self.__assign_tables_to_model_part(root_model_part, mat, table_dict)

        # Assign the properties to the model part's elements and conditions.
        for elem in model_part.Elements:
            current_number_props = root_model_part.NumberOfProperties()
            elem_props = root_model_part.GetProperties(current_number_props + 1, mesh_id)
            elem_props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law)

            elem.Properties = elem_props

            self.__assign_interpolated_properties(mat, table_dict, elem, elem_props)

        for cond in model_part.Conditions:
            current_number_props = root_model_part.NumberOfProperties()
            cond_props = root_model_part.GetProperties(current_number_props + 1, mesh_id)
            cond_props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law)

            cond.Properties = cond_props

            self.__assign_interpolated_properties(mat, table_dict, cond, cond_props)


    #### Private methods needed for the interpolation ####

    def __interpolation_is_required(self, data):
        """
        This function checks if at least one of the variables requires interpolation
        """
        for key, value in data["Material"]["Variables"].items():
            if self.__is_double_with_interpolation(value):
                return True
            elif self.__is_vector_with_interpolation(value):
                return True
            elif self.__is_matrix_with_interpolation(value):
                return True
            else:
                return False

    def __is_double_with_interpolation(self, parameter):
        """
        This function is the analogon to IsDouble(), but it checks if interpolation is required
        """
        return self.__has_interpolation_keyword(parameter)

    def __is_vector_with_interpolation(self, parameter):
        """
        This function is the analogon to IsVector(), but it checks if interpolation is required
        It does NOT throw if the Vector is not valid, type checking is done later
        """
        interpolation_required = False
        if not parameter.IsArray():
            return False

        if parameter.size() > 0:
            if not parameter[0].IsArray(): # then this could be a matrix
                return False

        for i in range(parameter.size()):
            if self.__has_interpolation_keyword(parameter[i]):
                interpolation_required = True
            else:
                if not parameter[i].IsNumber(): # this means that the vector is not valid
                    return False

        return interpolation_required

    def __is_matrix_with_interpolation(self, parameter):
        """
        This function is the analogon to IsMatrix(), but it checks if interpolation is required
        It does NOT throw if the Matrix is not valid, type checking is done later
        """
        interpolation_required = False
        if not parameter.IsArray():
            return False

        num_rows = parameter.size()
        if num_rows == 0: # parameter is an empty array/vector => "[]"
            return False

        for i in range(num_rows):
            row_i = parameter[i]
            if not row_i.IsArray():
                return False

            num_cols = row_i.size()
            if num_cols != parameter[0].size(): # num of cols is not consistent
                return False

            for j in range(num_cols):
                if self.__has_interpolation_keyword(row_i[j]):
                    interpolation_required = True
                else:
                    if not row_i[j].IsNumber(): # this means that the matrix is not valid
                        return False

        return interpolation_required


    def __has_interpolation_keyword(self, parameter):
        if parameter.IsSubParameter():
            if parameter.Has("@table"):
                return True
            else:
                raise ValueError("Variable Dict is not valid!") # TODO throw here or just return False?
        else:
            return False


    def __assign_tables_to_model_part(self, root_model_part, material_parameters, table_dict):
        for table_name in sorted(material_parameters["Tables"].keys()):
            if table_name in table_dict.keys():
                err_msg = 'Table names must be unique, trying to add: "' + table_name
                err_msg += '" which exists already!'
                raise NameError(err_msg)

            table_param = material_parameters["Tables"][table_name]

            input_variable = self._GetVariable(table_param["input_variable"].GetString())
            output_variable = self._GetVariable(table_param["output_variable"].GetString())

            self.__check_variable_type(input_variable, table_name)
            self.__check_variable_type(output_variable, table_name)

            if not table_param.Has("input_variable_location"):
                raise Exception("You need to specify a variable location for the interpolation!")

            input_variable_location = table_param["input_variable_location"].GetString()

            table = KratosMultiphysics.PiecewiseLinearTable()

            for i in range(table_param["data"].size()):
                table.AddRow(table_param["data"][i][0].GetDouble(), table_param["data"][i][1].GetDouble())

            table_id = root_model_part.NumberOfTables() + 1

            root_model_part.AddTable(table_id, table)

            table_info = {
                "table_id" : table_id,
                "table" : table,
                "input_variable" : input_variable,
                "output_variable" : output_variable,
                "input_variable_location" : input_variable_location }

            table_dict[table_name] = table_info

    def __assign_interpolated_properties(self, mat, table_dict, geom_entity, prop):
        # assign Tables (stored in the RootModelPart) to the properties
        for table_name, table_info in table_dict.items():
            input_variable = table_info["input_variable"]
            output_variable = table_info["output_variable"]
            table = table_info["table"]
            prop.SetTable(input_variable, output_variable, table)

        # assign the values to the properties and interpolate if needed
        for key, value in mat["Variables"].items():
            var = self._GetVariable(key)
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
            elif self.__is_double_with_interpolation(value):
                interpolated_double = self.__compute_interpolated_value(geom_entity, table_dict,
                                                                        key, value["@table"].GetString())
                prop.SetValue(var, interpolated_double)
            elif self.__is_vector_with_interpolation(value):
                pass
            elif self.__is_matrix_with_interpolation(value):
                pass
                # if value.Has("@table"):
                #     table_name = value["@table"].GetString()
                #     if table_name not in tables_for_interpolation.keys():
                #         raise KeyError("Table \"" + table_name + "\" needed for interpolation is not defined!")
                #     table_info = tables_for_interpolation[table_name]
                #     interpolated_value = self.__compute_interpolated_value(geom_entity, prop, table_info, var)
                #     prop.SetValue( var, interpolated_value )
                # else:
                #     raise ValueError("Variable keyword is not valid!")
            else:
                raise TypeError("Type of value is not available for " + key)








        err
        # tables_for_interpolation = {}
        # for key, table in mat["Tables"].items():
        #     input_var = self._GetItemFromModule(table["input_variable"].GetString())
        #     output_var = self._GetItemFromModule(table["output_variable"].GetString())
        #     self._SetTable(prop, table, input_var, output_var)
        #     table_name = key

        #     # save information abt the tables used for interpolation
        #     if table.Has("input_variable_location"):
        #         input_var_location = table["input_variable_location"].GetString()
        #         tables_for_interpolation[table_name] = [input_var,output_var,input_var_location]

        # Add / override the values of material parameters in the properties


    def __compute_interpolated_value(self, geom_entity, table_dict, variable_name, table_name):
        # Retrieve information needed for interpolation
        if table_name not in table_dict.keys():
            raise NameError('Table "' + table_name + '" not found')

        table_info = table_dict[table_name]

        input_variable = table_info["input_variable"]
        input_variable_location = table_info["input_variable_location"]
        table = table_info["table"]

        if input_variable_location == "geometric_entity":
            if geom_entity.Has(input_variable): # Values in Geom Entites are saved as Non-historical values (model_part_io.cpp)
                input_value = geom_entity.GetValue(input_variable)
            else:
                err_msg  = "Geometric Entity # " + str(geom_entity.Id)
                err_msg += " does not have " + input_variable.Name()
                raise ValueError(err_msg)
        elif input_variable_location == "nodes":
            nodes = geom_entity.GetNodes()
            input_value = 0
            for node in nodes:
                if node.SolutionStepsDataHas(input_variable): # Values in Nodes are saved as Historical values (model_part_io.cpp)
                    input_value += node.GetSolutionStepValue(input_variable)
                else:
                    err_msg  = "Node # " + str(node.Id)
                    err_msg += " does not have " + input_variable.Name()
                    raise ValueError(err_msg)
                input_value /= len(nodes)
        else:
            raise Exception("Type of input_variable_location \"" + input_variable_location + "\" is not valid!")

        interpolated_value = table.GetValue(input_value) # interpolate the value from the table

        return interpolated_value

    def __check_variable_type(self, variable, table_name):
        if(type(variable) != KratosMultiphysics.DoubleVariable and type(variable) != KratosMultiphysics.Array1DComponentVariable):
            err_msg = 'In table "' + table_name + '": Variable type of variable - '
            err_msg += variable.Name() + ' - is incorrect!\nMust be a scalar or a component'
            raise Exception(err_msg)

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

            self.__assign_interpolated_properties(mat, elem, elem_props)

        for cond in model_part.Conditions:
            current_number_props = model_part.NumberOfProperties()
            cond_props = model_part.GetProperties(current_number_props + 1, mesh_id)

            cond.Properties = cond_props

            self.__assign_interpolated_properties(mat, cond, cond_props)

    def __assign_interpolated_properties(self, mat, geom_entity, prop):
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
                    interpolated_value = self.__compute_interpolated_value(geom_entity, prop, table_info, var)
                    prop.SetValue( var, interpolated_value )
                else:
                    raise ValueError("Variable keyword is not valid!")
            else:
                raise ValueError("Type of value is not available")

    def __compute_interpolated_value(self, geom_entity, prop, table_info, dict_key):
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