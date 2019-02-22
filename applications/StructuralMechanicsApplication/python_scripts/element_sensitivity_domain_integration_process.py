from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.json_utilities import read_external_json

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ElementSensitivityDomainIntegrationProcess(model, settings["Parameters"])

class ElementSensitivityDomainIntegrationProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """
        This class integrates scalar element sensitivities (material and cross-section
        properties like CROSS_AREA or YOUNGS_MODULUS) within defined domains.
        The integration domains are defined by sub model parts of the sensitivity model part.
    """

    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the model contaning the model_parts
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "help"                             : "This class integrates element sensitivities within domains defined by sub model parts",
            "element_sensitivity_variables"    : [],
            "model_part_name"                  : "",
            "sensitivity_model_part_name"      : "",
            "sensitivity_sub_model_part_list"  : [],
            "time_frequency"                   : 1.00
        }""")

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_parameters)
        self.params = params
        self.model  = model

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        model_part_name = self.params["model_part_name"].GetString()
        sensitivity_model_part_name = self.params["sensitivity_model_part_name"].GetString()
        sensitivity_sub_model_part_names = [ self.params["sensitivity_sub_model_part_list"][i].GetString() for i in range( 0,self.params["sensitivity_sub_model_part_list"].size() ) ]

        # Get sensitivity model part
        if (sensitivity_model_part_name != ""):
            self.sensitivity_model_part = self.model[model_part_name].GetSubModelPart(sensitivity_model_part_name)
        else:
            self.sensitivity_model_part = self.model[model_part_name]

        # Get defined sub model parts of sensitivity model part as integration domains
        self.sensitivity_sub_model_parts = []
        if len(sensitivity_sub_model_part_names) is 0:
            self.sensitivity_sub_model_parts.append(self.sensitivity_model_part)
        else:
            for mp_name in sensitivity_sub_model_part_names:
                self.sensitivity_sub_model_parts.append(self.sensitivity_model_part.GetSubModelPart(mp_name))

        self.element_sensitivity_variables = self.__generate_variable_list_from_input(self.params["element_sensitivity_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()

        # Initialize counter
        self.time_counter = 0.0

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        dt = self.sensitivity_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt

        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            # loop over integration domains
            for sub_mp_i in self.sensitivity_sub_model_parts:
                # loop over sensitivty variables for which integration should performed
                for variable_i in self.element_sensitivity_variables:
                    value = 0.0
                    # loop over elements in sensitivity sub model parts
                    if next(sub_mp_i.Elements.__iter__() ).Has(variable_i):
                        for elem_i in sub_mp_i.Elements:
                            value += elem_i.GetValue(variable_i)
                        for elem_i in sub_mp_i.Elements:
                            elem_i.SetValue(variable_i, value)
                    else:
                        raise Exception(variable_i.Name() + " is not available for domain integration!")

    def __generate_variable_list_from_input(self,param):
        """ Parse a list of variables from input.

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        variable_list = []
        for i in range(0, param.size()):
            variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(param[i].GetString())
            if (variable_type == "Double"):
                variable_list.append(KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ))
            else:
                raise Exception("sensitivity domain integration is only available for variables of data type 'Double' but " + param[i].GetString() + " is of type '" + variable_type + "'.")

        return variable_list
