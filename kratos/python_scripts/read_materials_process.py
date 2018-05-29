from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
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

        See __AssignPropertyBlock for detail on how properties are imported.
        """
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "materials_filename" : "please specify the file to be opened"
        }
        """)

        settings.ValidateAndAssignDefaults(default_settings)
        self.Model = Model

        KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Started")

        with open(settings["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        for i in range(materials["properties"].size()):
            data = materials["properties"][i]

            if InterpolationIsRequired(data):
                self.__AssignPropertyBlockInterpolated(data)
            else:
                self.__AssignPropertyBlock(data)

        KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Finished")

    def __AssignPropertyBlock(self, data):
        """Set constitutive law and material properties and assign to elements and conditions.

        Arguments:
        data -- a dictionary or json object defining properties for a model part.

        Example:
        data = {
            "model_part_name" : "Plate",
            "properties_id" : 1,
            "Material" : {
                "constitutive_law" : {
                    "name" : "LinearElasticPlaneStress2DLaw"
                },
                "Variables" : {
                    "YOUNG_MODULUS" : 200e9,
                    "POISSON_RATIO" : 0.3,
                    "RESIDUAL_VECTOR" : [1.5,0.3,-2.58],
                    "LOCAL_INERTIA_TENSOR" : [[0.27,0.0],[0.0,0.27]]
                },
                "Tables" : {
                    "Table1" : {
                        "input_variable" : "TEMPERATURE",
                        "output_variable" : "YOUNG_MODULUS",
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
        property_id = data["properties_id"].GetInt()
        mesh_id = 0
        prop = model_part.GetProperties(property_id, mesh_id)

        if len(data["Material"]["Variables"].keys()) > 0 and prop.HasVariables():
            KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Property", str(property_id), "already has variables." )
        if len(data["Material"]["Tables"].keys()) > 0 and prop.HasTables():
            KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Property", str(property_id), "already has tables." )

        # Assign the properties to the model part's elements and conditions.
        for elem in model_part.Elements:
            elem.Properties = prop

        for cond in model_part.Conditions:
            cond.Properties = prop

        mat = data["Material"]

        # Set the CONSTITUTIVE_LAW for the current properties.
        constitutive_law = GetConstitutiveLaw( mat["constitutive_law"]["name"].GetString() )

        prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law.Clone())

        # Add / override the values of material parameters in the properties
        for var_name, value in mat["Variables"].items():
            var = GetVariable(var_name)
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
                raise TypeError("Type of value is not available for " + var_name)

        # Add / override tables in the properties
        for table_name, table_data in mat["Tables"].items():
            input_variable = GetVariable(table_data["input_variable"].GetString())
            output_variable = GetVariable(table_data["output_variable"].GetString())

            CheckTableVariableType(input_variable, table_name)
            CheckTableVariableType(output_variable, table_name)

            table = GetTable(table_data)

            prop.SetTable(input_variable, output_variable, table)

    def __AssignPropertyBlockInterpolated(self, data):
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
        constitutive_law = GetConstitutiveLaw( mat["constitutive_law"]["name"].GetString() )

        # Set the tables to the ModelPart
        table_dict = {}
        AssignTablesToModelPart(root_model_part, mat, table_dict)

        # Assign the properties to the model part's elements and conditions.
        for elem in model_part.Elements:
            current_number_props = root_model_part.NumberOfProperties()
            elem_props = root_model_part.GetProperties(current_number_props + 1, mesh_id)
            elem_props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law.Clone())

            elem.Properties = elem_props

            AssignInterpolatedProperties(mat, table_dict, elem, elem_props)

        for cond in model_part.Conditions:
            current_number_props = root_model_part.NumberOfProperties()
            cond_props = root_model_part.GetProperties(current_number_props + 1, mesh_id)
            cond_props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law.Clone())

            cond.Properties = cond_props

            AssignInterpolatedProperties(mat, table_dict, cond, cond_props)

def GetAttribute(my_string, function_pointer, attribute_type):
    """Return the python object named by the string argument.

    To be used with functions from KratosGlobals

    Examples:
    variable = GetAttribute("DISPLACEMENT",
                            KratosMultiphysics.KratosGlobals.GetVariable,
                            "Variable")

    constitutive_law = GetAttribute("LinearElastic3DLaw",
                                    KratosMultiphysics.KratosGlobals.GetConstitutiveLaw,
                                    "Constitutive Law")
    """
    splitted = my_string.split(".")

    if len(splitted) == 0:
        raise Exception("Something wrong. Trying to split the string " + my_string)
    if len(splitted) > 3:
        raise Exception("Something wrong. String " + my_string + " has too many arguments")

    attribute_name = splitted[-1]

    if len(splitted) == 2 or len(splitted) == 3:
        warning_msg =  "Ignoring \"" +  my_string.rsplit(".",1)[0]
        warning_msg += "\" for " + attribute_type +" \"" + attribute_name + "\""
        KratosMultiphysics.Logger.PrintInfo("Warning in reading materials", warning_msg)

    return function_pointer(attribute_name) # This also checks if the application has been imported

def GetVariable(my_string):
    """Return the python object of a Variable named by the string argument.

    Examples:
    recommended usage:
    variable = GetVariable("VELOCITY")
    deprecated:
    variable = GetVariable("KratosMultiphysics.VELOCITY")
    variable = GetVariable("SUBSCALE_PRESSURE")
    variable = GetVariable("FluidDynamicsApplication.SUBSCALE_PRESSURE")
    variable = GetVariable("KratosMultiphysics.FluidDynamicsApplication.SUBSCALE_PRESSURE")
    """
    return GetAttribute(my_string, KratosMultiphysics.KratosGlobals.GetVariable, "Variable")

def GetConstitutiveLaw(my_string):
    """Return the python object of a Constitutive Law named by the string argument.

    Example:
    recommended usage:
    constitutive_law = GetConstitutiveLaw('LinearElastic3DLaw')
    deprecated:
    constitutive_law = GetConstitutiveLaw('KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw')
    constitutive_law = GetConstitutiveLaw('StructuralMechanicsApplication.LinearElastic3DLaw')

    model_part.GetProperties(prop_id).SetValue(CONSTITUTIVE_LAW, constitutive_law)
    """
    return GetAttribute(my_string, KratosMultiphysics.KratosGlobals.GetConstitutiveLaw, "Constitutive Law")

    #### Methods needed for the interpolation ####

def InterpolationIsRequired(data):
    """
    This function checks if at least one of the variables requires interpolation
    """
    for var_name, value in data["Material"]["Variables"].items():
        if IsDoubleWithInterpolation(value):
            return True
        elif IsVectorWithInterpolation(value):
            return True
        elif IsMatrixWithInterpolation(value):
            return True

    return False

def IsDoubleWithInterpolation(parameter):
    """
    This function is the analogon to IsDouble(), but it checks if interpolation is required
    """
    return HasInterpolationKeyword(parameter)

def IsVectorWithInterpolation(parameter):
    """
    This function is the analogon to IsVector(), but it checks if interpolation is required
    It does NOT throw if the Vector is not valid, type checking is done later
    """
    interpolation_required = False
    if not parameter.IsArray():
        return False

    if parameter.size() > 0:
        if parameter[0].IsArray(): # then this could be a matrix
            return False

    for i in range(parameter.size()):
        if HasInterpolationKeyword(parameter[i]):
            interpolation_required = True
        else:
            if not parameter[i].IsNumber(): # this means that the vector is not valid
                return False

    return interpolation_required

def IsMatrixWithInterpolation(parameter):
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
            if HasInterpolationKeyword(row_i[j]):
                interpolation_required = True
            else:
                if not row_i[j].IsNumber(): # this means that the matrix is not valid
                    return False

    return interpolation_required

def HasInterpolationKeyword(parameter):
    """
    This function checks if the parameter has the keyword for Interpolation
    """
    if parameter.IsSubParameter():
        if parameter.Has("@table"):
            return True
        else:
            return False
    else:
        return False

def SizeInterpolatedVector(vector_parameter):
    """
    This function returns the size of a vector that requires Interpolation
    """
    if vector_parameter.IsVector() or IsVectorWithInterpolation(vector_parameter):
        return vector_parameter.size()
    else:
        raise TypeError("Object is not a Vector!")

def SizeInterpolatedMatrix(matrix_parameter):
    """
    This function returns the size of a matrix that requires Interpolation
    """
    if matrix_parameter.IsMatrix() or IsMatrixWithInterpolation(matrix_parameter):
        # Existance of these values is assured through the checks above
        size1 = matrix_parameter.size()
        size2 = matrix_parameter[0].size()
        return size1, size2
    else:
        raise TypeError("Object is not a Matrix!")

def AssignTablesToModelPart(root_model_part, material_parameters, table_dict):
    """
    This function reads the tables and assigns them to the modelpart
    """
    for table_name in sorted(material_parameters["Tables"].keys()):
        if table_name in table_dict.keys():
            err_msg = 'Table names must be unique, trying to add: "' + table_name
            err_msg += '" which exists already!'
            raise NameError(err_msg)

        table_data = material_parameters["Tables"][table_name]

        input_variable = GetVariable(table_data["input_variable"].GetString())
        output_variable = GetVariable(table_data["output_variable"].GetString())

        CheckTableVariableType(input_variable, table_name)
        CheckTableVariableType(output_variable, table_name)

        if not table_data.Has("input_variable_location"):
            raise Exception("You need to specify a variable location for the interpolation!")

        input_variable_location = table_data["input_variable_location"].GetString()

        table = GetTable(table_data)

        table_id = root_model_part.NumberOfTables() + 1

        root_model_part.AddTable(table_id, table)

        table_info = {
            "table_id" : table_id,
            "table" : table,
            "input_variable" : input_variable,
            "output_variable" : output_variable,
            "input_variable_location" : input_variable_location }

        table_dict[table_name] = table_info

def AssignInterpolatedProperties(mat, table_dict, geom_entity, prop):
    """
    This function assigns the interpolated quantities to the properties
    Note that the table is stored in the modelpart and passed to the property
    """
    # assign Tables (stored in the RootModelPart) to the properties
    for table_name, table_info in table_dict.items():
        input_variable = table_info["input_variable"]
        output_variable = table_info["output_variable"]
        table = table_info["table"]
        prop.SetTable(input_variable, output_variable, table)

    # assign the values to the properties and interpolate if needed
    for var_name, value in mat["Variables"].items():
        var = GetVariable(var_name)
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
        elif IsDoubleWithInterpolation(value):
            interpolated_double = ComputeInterpolatedValue(geom_entity, table_dict,
                                                            var_name, value)
            prop.SetValue(var, interpolated_double)
        elif IsVectorWithInterpolation(value):
            interpolated_vector = ComputeInterpolatedVector(geom_entity, table_dict,
                                                            var_name, value)
            prop.SetValue(var, interpolated_vector)
        elif IsMatrixWithInterpolation(value):
            interpolated_matrix = ComputeInterpolatedMatrix(geom_entity, table_dict,
                                                            var_name, value)
            prop.SetValue(var, interpolated_matrix)
        else:
            raise TypeError("Type of value is not available for " + var_name)

def ComputeInterpolatedValue(geom_entity, table_dict, variable_name, value):
    """
    This function interpolates the value
    the input value can be on "nodes" or on "geom_entity"
    """
    # Retrieve information needed for interpolation
    table_name = value["@table"].GetString()

    if table_name not in table_dict.keys():
        raise NameError('Table "' + table_name + '" not found')

    table_info = table_dict[table_name]

    input_variable = table_info["input_variable"]
    input_variable_location = table_info["input_variable_location"]
    table = table_info["table"]

    if input_variable_location == "geom_entity":
        if geom_entity.Has(input_variable): # Values in Geom Entites are saved as Non-historical values (model_part_io.cpp)
            input_value = geom_entity.GetValue(input_variable) # This is a double, since Tables exist only with doubles!
        else:
            err_msg  = "Geometric Entity # " + str(geom_entity.Id)
            err_msg += " does not have " + input_variable.Name()
            raise ValueError(err_msg)
    elif input_variable_location == "nodes":
        nodes = geom_entity.GetNodes()
        input_value = 0.0
        for node in nodes:
            if node.SolutionStepsDataHas(input_variable): # Values in Nodes are saved as Historical values (model_part_io.cpp)
                input_value += node.GetSolutionStepValue(input_variable) # This is a double, since Tables exist only with doubles!
            else:
                err_msg  = "Node # " + str(node.Id)
                err_msg += " does not have " + input_variable.Name()
                raise ValueError(err_msg)
        input_value /= len(nodes)
    else:
        raise Exception('Type of input_variable_location "' + input_variable_location + '" is not valid!')

    interpolated_value = table.GetValue(input_value) # interpolate the value from the table

    return interpolated_value

def ComputeInterpolatedVector(geom_entity, table_dict, variable_name, vector_parameter):
    """
    This function computes a vector where some values require interpolation
    """
    size = SizeInterpolatedVector(vector_parameter)
    interpolated_vector = KratosMultiphysics.Vector(size)

    for i in range(size):
        sub_param = vector_parameter[i]
        if sub_param.IsDouble():
            interpolated_vector[i] = sub_param.GetDouble()
        elif sub_param.IsInt(): # needed?
            interpolated_vector[i] = sub_param.GetInt()
        elif HasInterpolationKeyword(sub_param):
            interpolated_vector[i] = ComputeInterpolatedValue(geom_entity, table_dict,
                                                              variable_name, sub_param)
        else:
            raise TypeError("Wrong Type of Value")

    return interpolated_vector

def ComputeInterpolatedMatrix(geom_entity, table_dict, variable_name, matrix_parameter):
    """
    This function computes a matrix where some values require interpolation
    """
    size1, size2 = SizeInterpolatedMatrix(matrix_parameter)
    interpolated_matrix = KratosMultiphysics.Matrix(size1, size2)

    for i in range(size1):
        for j in range(size2):
            sub_param = matrix_parameter[i][j]
            if sub_param.IsDouble():
                interpolated_matrix[i,j] = sub_param.GetDouble()
            elif sub_param.IsInt(): # needed?
                interpolated_matrix[i,j] = sub_param.GetInt()
            elif HasInterpolationKeyword(sub_param):
                interpolated_matrix[i,j] = ComputeInterpolatedValue(geom_entity, table_dict,
                                                                    variable_name, sub_param)
            else:
                raise TypeError("Wrong Type of Value")

    return interpolated_matrix

def CheckTableVariableType(variable, table_name):
    """
    This function checks if the variable type is suitable
    """
    var_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable.Name())
    if(var_type != "Double" and var_type != "Component"):
        err_msg = 'In table "' + table_name + '": Variable type of variable - '
        err_msg += variable.Name() + ' - is incorrect!\nMust be a scalar or a component'
        raise TypeError(err_msg)

def GetTable(table_data):
    """
    This function reads a table from a Kratos::Parameters Object
    TODO move this function to Kratos::Parameters
    """
    table = KratosMultiphysics.PiecewiseLinearTable()

    for i in range(table_data["data"].size()):
        table.AddRow(table_data["data"][i][0].GetDouble(), table_data["data"][i][1].GetDouble())

    return table