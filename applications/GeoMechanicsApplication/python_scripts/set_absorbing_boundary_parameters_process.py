from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetAbsorbingBoundaryParametersProcess(Model, settings["Parameters"])


class SetAbsorbingBoundaryParametersProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("absorbing_factors", settings["absorbing_factors"])
        params.AddValue("virtual_thickness", settings["virtual_thickness"])

        self.process = KratosGeo.SetAbsorbingBoundaryParametersProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()
