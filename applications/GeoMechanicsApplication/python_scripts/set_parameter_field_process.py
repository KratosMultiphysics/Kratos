import json
import importlib

import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    # return KratosGeo.SetParameterFieldProcess(Model[settings["Parameters"]["model_part_name"].GetString()], settings["Parameters"])
    return SetParameterFieldProcess(Model, settings["Parameters"])
#
class SetParameterFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]


        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.params.AddValue("variable_name", settings["variable_name"])
        self.params.AddValue("func_type", settings["func_type"])
        self.params.AddValue("function", settings["function"])
        self.params.AddValue("dataset", settings["dataset"])
        self.process = KratosGeo.SetParameterFieldProcess(self.model_part, self.params)
#
    def ExecuteInitialize(self):

        if self.params["func_type"].GetString() == "python":
            elements = [element for element in self.model_part.Elements]

            input_dict = {}
            return_dict = {}
            i=0

            for element in self.model_part.Elements:
            #for i in range(len(elements)):
                center = element.GetGeometry().Center()
                coords = [center[0], center[1], center[2]]
                value = element.Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)


                #element_dict["id"] = elements[i].Id
                #element_dict["variable"] = "YOUNG_MODULUS"
                #element_dict["value"] = value
                #element_dict["coordinates"] = coords
                return_dict[str(elements[i].Id)] = value
                input_dict[str(element.Id)] ={"value": value, "coordinates": coords}
                # input_dict[str(element.Id)]["value"] = value
                # input_dict[str(element.Id)]["coordinates"] = coords
                #return_dict[str(element.Id)] = value

                #i+=1

                #properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, 500 *i)

                #elements[i].Properties = properties
                #E2 = elements[i].Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
                #a=1+1

            custom_module = importlib.import_module("."+self.params["function"].GetString(),KratosGeo.__name__)
            custom_module.run(input_dict, return_dict)

            self.params["dataset"].SetString(json.dumps(return_dict))

        else:
            self.process.ExecuteInitialize()

            # self.params.AddValue("time", KratosMultiphysics.KratosGlobals.GetVariable("TIME"))
            # KratosGeo.SetParameterFieldProcess(self.model_part, self.params)

        self.process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()