import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.CableNetApplication as CableNetApplication


from numpy import polyfit

from KratosMultiphysics import Logger



#/**
# * @class EmpiricalSpringElementProcess
# *
# * @brief This process creates a spring element w.r.t. to given displacement/load data points
# *
# * @author Klaus B Sautter
# */



def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EmpiricalSpringElementProcess(Model, settings["Parameters"])



class EmpiricalSpringElementProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"           : "example_part",
            "computing_model_part_name" : "Structure",
            "node_ids"                  : [1,2],
            "element_id"                : 1,
            "property_id"               : 1,
            "displacement_data"         : [0.0,1.0,2.0,3.0],
            "force_data"                : [0.0,1.0,2.0,3.0],
            "polynomial_order"          : 3
        }
        """)
        default_settings.ValidateAndAssignDefaults(settings)
        self.function_fitted = polyfit(settings["displacement_data"].GetVector(),settings["force_data"].GetVector(),settings["polynomial_order"].GetInt())


        # The computing model part
        self.computing_model_part = Model[settings["computing_model_part_name"].GetString()]
        self.custom_model_part     = Model[settings["model_part_name"].GetString()]

        self.empirical_spring_element_process = CableNetApplication.EmpiricalSpringElementProcess(self.custom_model_part, settings, self.function_fitted)


    def ExecuteInitialize(self):
        self.empirical_spring_element_process.ExecuteInitialize()

        ## add new element in the computing MP
        for element_i in self.custom_model_part.Elements:
            self.computing_model_part.AddElement(element_i, 0)
        Logger.PrintInfo("Initialized","EmpiricalSpringElementProcess")

