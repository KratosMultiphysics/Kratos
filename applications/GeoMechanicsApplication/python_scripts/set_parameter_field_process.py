import json
import importlib

import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
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

            values = []
            all_coordinates = []

            for element in self.model_part.Elements:
                center = element.GetGeometry().Center()
                coords = [center[0], center[1], center[2]]
                value = element.Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)

                values.append(value)
                all_coordinates.append(coords)

                return_dict[str(elements[i].Id)] = value
                #input_dict[str(element.Id)] = {"value": value, "coordinates": coords}

            input_dict["values"] = values
            input_dict["coordinates"] = all_coordinates

            custom_module = importlib.import_module("." + self.params["function"].GetString(), KratosGeo.__name__ +
                                                    ".user_defined_scripts")
            custom_module.run(input_dict, return_dict)

            self.params["dataset"].SetString(json.dumps(return_dict))

        self.process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()