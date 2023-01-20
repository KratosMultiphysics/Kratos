
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return KratosGeo.SetParameterFieldProcess(Model[settings["Parameters"]["model_part_name"].GetString()], settings["Parameters"])

#
# class SetParameterFieldProcess(KratosMultiphysics.Process):
#     def __init__(self, Model, settings ):
#         KratosMultiphysics.Process.__init__(self)
#
#         model_part = Model[settings["model_part_name"].GetString()]
#
#         params = KratosMultiphysics.Parameters("{}")
#         params.AddValue("model_part_name",settings["model_part_name"])
#         params.AddValue("variable_name", settings["variable_name"])
#         params.AddValue("func_type", settings["func_type"])
#         params.AddValue("function", settings["function"])
#         KratosGeo.SetParameterFieldProcess(model_part, params)
#
#     def ExecuteInitialize(self):
#         self.process.ExecuteInitialize()