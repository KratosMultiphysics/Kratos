from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import sys

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ExternalSensitivityProcess(model, settings["Parameters"])

def CheckAvailabilityOfElementSensitivities(variable, model_part):
    """ Check if element sensitivities w.r.t. given variable are available.
    Keyword arguments:
    variable -- traced variable within sensitivity analysis.
    model_part -- sub model part of the sensitivity model part.
    """

    if sys.version_info[0] >= 3: # python3 syntax
        return model_part.Elements.__iter__().__next__().Has(variable)
    else: # python2 syntax
        return model_part.Elements.__iter__().next().Has(variable)

def CheckAvailabilityOfNodalSensitivities(variable, model_part):
    """ Check if element sensitivities w.r.t. given variable are available.
    Keyword arguments:
    variable -- traced variable within sensitivity analysis.
    model_part -- sub model part of the sensitivity model part.
    """

    if sys.version_info[0] >= 3: # python3 syntax
        return model_part.Nodes.__iter__().__next__().Has(variable)
    else: # python2 syntax
        return model_part.Nodes.__iter__().next().Has(variable)

def CheckAvailabilityOfConditionSensitivities(variable, model_part):
    """ Check if element sensitivities w.r.t. given variable are available.
    Keyword arguments:
    variable -- traced variable within sensitivity analysis.
    model_part -- sub model part of the sensitivity model part.
    """

    if sys.version_info[0] >= 3: # python3 syntax
        return model_part.Conditions.__iter__().__next__().Has(variable)
    else: # python2 syntax
        return model_part.Conditions.__iter__().next().Has(variable)


class ExternalSensitivityProcess(KratosMultiphysics.Process):
    """
        This class integrates scalar element sensitivities (material and cross-section
        properties like CROSS_AREA or YOUNGS_MODULUS) within defined domains.
        The integration domains are defined by sub model parts of the sensitivity model part.
    """

    def __init__(self, model, parameter):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the model contaning the model_parts
        parameter -- Kratos parameters containing process settings.
        """
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "help"                             : "This process initializes the non-historical data base for external computed sensitivities",
            "model_part_name"                             : "",
            "nodal_data_value_sensitivity_variables"      : [],
            "element_data_value_sensitivity_variables"    : [],
            "condition_data_value_sensitivity_variables"  : []
        }""")

        ## Overwrite the default settings with user-provided parameters
        parameter.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[parameter["model_part_name"].GetString()]

        self.nodal_sensitivity_variables = self.__GenerateVariableListFromInput(parameter["nodal_data_value_sensitivity_variables"])
        self.element_sensitivity_variables = self.__GenerateVariableListFromInput(parameter["element_data_value_sensitivity_variables"])
        self.condition_sensitivity_variables = self.__GenerateVariableListFromInput(parameter["condition_data_value_sensitivity_variables"])


    def ExecuteInitialize(self):
        # initialize all non-hinstorical external sensitivity variables
        for variable_i in self.nodal_sensitivity_variables:
            if not CheckAvailabilityOfNodalSensitivities(variable_i, self.model_part):
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(variable_i, self.model_part.Nodes)
        for variable_i in self.element_sensitivity_variables:
            if not CheckAvailabilityOfElementSensitivities(variable_i, self.model_part):
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(variable_i, self.model_part.Elements)
        for variable_i in self.condition_sensitivity_variables:
            if not CheckAvailabilityOfConditionSensitivities(variable_i, self.model_part):
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(variable_i, self.model_part.Conditions)


    def __GenerateVariableListFromInput(self, parameter):
        """ Parse a list of variables from input.
        Keyword arguments:
        self -- It signifies an instance of a class.
        parameter -- Kratos parameters containing process settings.
        """

        if not parameter.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        variable_list = []
        for i in range(0, parameter.size()):
            variable_list.append(KratosMultiphysics.KratosGlobals.GetVariable( parameter[i].GetString() ))

        return variable_list







