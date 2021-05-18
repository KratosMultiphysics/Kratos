# Importing the Kratos Library
import KratosMultiphysics
import sys

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ElementSensitivityDomainIntegrationProcess(model, settings["Parameters"])

def CheckAvailabilityOfSensitivities(variable, model_part):
    """ Check if element sensitivities w.r.t. given variable are available.
    Keyword arguments:
    variable -- traced variable within sensitivity analysis.
    model_part -- sub model part of the sensitivity model part.
    """
    return model_part.Elements.__iter__().__next__().Has(variable)

class ElementSensitivityDomainIntegrationProcess(KratosMultiphysics.Process):
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
            "help"                             : "This class integrates element sensitivities within domains defined by sub model parts",
            "element_sensitivity_variables"    : [],
            "model_part_name"                  : "",
            "sensitivity_model_part_name"      : "",
            "sensitivity_sub_model_part_list"  : []
        }""")

        ## Overwrite the default settings with user-provided parameters
        parameter.ValidateAndAssignDefaults(default_parameters)

        model_part_name = parameter["model_part_name"].GetString()
        sensitivity_model_part_name = parameter["sensitivity_model_part_name"].GetString()
        sensitivity_sub_model_part_names = [ parameter["sensitivity_sub_model_part_list"][i].GetString() for i in range( 0, parameter["sensitivity_sub_model_part_list"].size() ) ]

        # Get sensitivity model part
        if (sensitivity_model_part_name != ""):
            self.sensitivity_model_part = model[model_part_name].GetSubModelPart(sensitivity_model_part_name)
        else:
            self.sensitivity_model_part = model[model_part_name]

        # Get defined sub model parts of sensitivity model part as integration domains
        self.sensitivity_sub_model_parts = []
        if len(sensitivity_sub_model_part_names) is 0:
            self.sensitivity_sub_model_parts.append(self.sensitivity_model_part)
        else:
            for mp_name in sensitivity_sub_model_part_names:
                self.sensitivity_sub_model_parts.append(self.sensitivity_model_part.GetSubModelPart(mp_name))

        self.element_sensitivity_variables = self.__GenerateVariableListFromInput(parameter["element_sensitivity_variables"])


    def Check(self):
        """ This method is executed at the begining to verify that the input is correct.
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Check integration domains
        for sub_mp_i in self.sensitivity_sub_model_parts:
            if sub_mp_i.NumberOfElements() < 1:
                raise Exception("sensitivity sub model part has no elements!")
        return 0


    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # loop over sensitivty variables for which integration should performed
        for variable_i in self.element_sensitivity_variables:
            if CheckAvailabilityOfSensitivities(variable_i, self.sensitivity_model_part):
                # loop over integration domains
                for sub_mp_i in self.sensitivity_sub_model_parts:
                    domain_sensitivity = KratosMultiphysics.VariableUtils().SumElementScalarVariable(variable_i, sub_mp_i)
                    KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(variable_i, domain_sensitivity, sub_mp_i.Elements)
            else:
                raise Exception(variable_i.Name() + " is not available for domain integration!")


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
            variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(parameter[i].GetString())
            if (variable_type == "Double"):
                variable_list.append(KratosMultiphysics.KratosGlobals.GetVariable( parameter[i].GetString() ))
            else:
                raise Exception("sensitivity domain integration is only available for variables of data type 'Double' but " + parameter[i].GetString() + " is of type '" + variable_type + "'.")

        return variable_list






