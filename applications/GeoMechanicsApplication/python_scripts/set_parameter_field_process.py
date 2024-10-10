import json
import importlib

import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo


def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return SetParameterFieldProcess(Model, settings["Parameters"])


class SetParameterFieldProcess(KratosMultiphysics.Process):
    """
    Sets parameter field process. This process has 3 option to generate a custom parameter field:

    | option 1, func_type = 'input': with this option, a parameter field is generated based on a function which is
    directly written in the projectparameters.json at the 'function' parameter. This function can depend on
    the x, y and z coordinate.

    | option 2, func_type = 'python': with this option, a parameter field is generated, using a user defined python
    script. This python script has to be inherited from the 'ParameterFieldBase' class. Which is located at:
    'GeoMechanicsApplication->python_scripts->user_defined_scripts->user_defined_parameter_field_base.py'
    the name of the script (without '.py') should be filled in at the 'function' parameter in the projectparameters.json

    | option 3, func_type = 'json_file': with this option, a parameter field can be directly read from a json dictionary.
    This dictionary has to contain the 'values' key, which is a 1D list of all the field values. The list has to have
    the same size as the elements within the model part, and need to be accordingly sorted. The filename should be
    filled in at the 'dataset' parameter, within the projectparameters.json
    """

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name", settings["model_part_name"])
        self.params.AddValue("variable_name", settings["variable_name"])
        self.params.AddValue("func_type", settings["func_type"])
        self.params.AddValue("function", settings["function"])
        self.params.AddValue("dataset", settings["dataset"])
        is_variable_vector_type = isinstance(self.GetVariableBasedOnString(), KratosMultiphysics.VectorVariable)
        if "input" in settings["func_type"].GetString() and is_variable_vector_type:
            self.params.AddValue("vector_variable_indices", settings["vector_variable_indices"])
        if "json_file" in settings["func_type"].GetString():
            self.params.AddValue("dataset_file_name", settings["dataset_file_name"])
        self.process = KratosGeo.SetParameterFieldProcess(self.model_part, self.params)

    def GetVariableBasedOnString(self):
        """
        This function returns the variable based on the variable name string.


        Returns
        -------
        variable : KratosMultiphysics.Variable

        """

        # Get variable object
        imported_modules = [KratosGeo, KratosMultiphysics]

        for kratos_module in imported_modules:
            if hasattr(kratos_module, self.params["variable_name"].GetString()):
                variable = getattr(kratos_module, self.params["variable_name"].GetString())
                return variable

        raise AttributeError(f'The variable: {self.params["variable_name"].GetString()} is not present within '
                             f'the imported modules')


    def ExecuteInitialize(self):
        """
        Initializes the process. Within the python part of 'ExecuteInitialize', the parameter field with func type
        'python' is generate. The cpp part sets the fields on the elements.

        Returns
        -------

        """

        # if the type of parameter field is a user defined python script:
        if self.params["func_type"].GetString() == "python":

            self.params["func_type"].SetString("json_string")

            # initialise input and output
            input_dict = {}
            return_dict = {}

            values = []
            all_coordinates = []

            variable = self.GetVariableBasedOnString()

            for element in self.model_part.Elements:

                # calculate center coordinates of the element
                center = element.GetGeometry().Center()
                coords = [center[0], center[1], center[2]]

                # Get value of parameter at current element
                value = element.Properties.GetValue(variable)

                values.append(value)
                all_coordinates.append(coords)

            # add variable values and coordinates to input dictionary
            input_dict["values"] = values
            input_dict["coordinates"] = all_coordinates

            # import user defined module
            custom_module = importlib.import_module("." + self.params["function"].GetString(),
                                                                KratosGeo.__name__ + ".user_defined_scripts")

            CustomParameterField = getattr(custom_module, 'ParameterField')
            custom_class = CustomParameterField()

            # validate and generate custom defined parameter field
            custom_class.validate_input(input_dict, return_dict)
            custom_class.generate_field()
            custom_class.validate_output()

            self.params["dataset"].SetString(json.dumps(return_dict))

        # run cpp part of ExecuteInitialize
        self.process.ExecuteInitialize()
