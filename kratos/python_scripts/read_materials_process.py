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

        See _AssignPropertyBlock for detail on how properties are imported.
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
            self._AssignPropertyBlock(materials["properties"][i])

        KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Finished")

    def _get_attribute(self, my_string, function_pointer, attribute_type):
        """Return the python object named by the string argument.

        To be used with functions from KratosGlobals

        Examples:
        variable = self._get_attribute("DISPLACEMENT",
                                       KratosMultiphysics.KratosGlobals.GetVariable,
                                       "Variable")

        constitutive_law = self._get_attribute("LinearElastic3DLaw",
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


    def _GetVariable(self, my_string):
        """Return the python object of a Variable named by the string argument.

        Examples:
        recommended usage:
        variable = self._GetVariable("VELOCITY")
        deprecated:
        variable = self._GetVariable("KratosMultiphysics.VELOCITY")
        variable = self._GetVariable("SUBSCALE_PRESSURE")
        variable = self._GetVariable("FluidDynamicsApplication.SUBSCALE_PRESSURE")
        variable = self._GetVariable("KratosMultiphysics.FluidDynamicsApplication.SUBSCALE_PRESSURE")
        """
        return self._get_attribute(my_string, KratosMultiphysics.KratosGlobals.GetVariable, "Variable")

    def _GetConstitutiveLaw(self, param):
        """Return the python object of a Constitutive Law named by the string argument.

        Example:
        recommended usage:
        constitutive_law = self._GetConstitutiveLaw('LinearElastic3DLaw')
        deprecated:
        constitutive_law = self._GetConstitutiveLaw('KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw')
        constitutive_law = self._GetConstitutiveLaw('StructuralMechanicsApplication.LinearElastic3DLaw')

        model_part.GetProperties(prop_id).SetValue(CONSTITUTIVE_LAW, constitutive_law)
        """
        my_string = param["name"].GetString()
        splitted = my_string.split(".")

        if len(splitted) == 0:
            raise Exception("Something wrong. Trying to split the string " + my_string)
        if len(splitted) > 3:
            raise Exception("Something wrong. String " + my_string + " has too many arguments")

        cl_name = splitted[-1]
        param["name"].SetString(cl_name)
        cl = self._get_attribute(cl_name, KratosMultiphysics.KratosGlobals.GetConstitutiveLaw, "Constitutive Law")
        return cl.Create(param)

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
        if (mat.Has("constitutive_law")):
            constitutive_law = self._GetConstitutiveLaw( mat["constitutive_law"] )

            prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law.Clone())
        else:
            KratosMultiphysics.Logger.PrintWarning("::[Reading materials process]:: ", "Not constitutive law defined for material ID: ", property_id)

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
