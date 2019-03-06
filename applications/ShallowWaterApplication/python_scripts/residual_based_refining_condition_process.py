from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ResidualBasedRefiningConditionProcess(Model, settings["Parameters"])

class ResidualBasedRefiningConditionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"           : "model_part",
            "error_variable"            : "RESIDUAL_NORM",
            "variable_threshold"        : 1e-3,
            "only_refine_wet_domain"    : true
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.error_variable = getattr(KratosMultiphysics, settings["error_variable"].GetString())
        self.variable_threshold = settings["variable_threshold"].GetDouble()
        self.only_refine_wet_domain = settings["only_refine_wet_domain"].GetBool()

    def Execute(self):
        # TODO: move to a c++ process
        for node in self.model_part.Nodes:
            node.Set(KratosMultiphysics.TO_REFINE, False)
        # Set the nodes or the elements which are to refine
        for elem in self.model_part.Elements:
            residual = elem.GetValue(KratosMultiphysics.RESIDUAL_NORM)
            if residual > self.variable_threshold:
                active_element = True
                if self.only_refine_wet_domain:
                    element_is_wet = False
                    for node in elem.GetNodes():
                        if node.GetSolutionStepValue(Shallow.HEIGHT) > 0.0:
                            element_is_wet = True
                    active_element = element_is_wet
                if active_element:
                    for node in elem.GetNodes():
                        node.Set(KratosMultiphysics.TO_REFINE, True)
