import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarConstraintTableProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyScalarConstraintTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.params.AddValue("variable_name",settings["variable_name"])
        if settings.Has("is_fixed"):
            self.params.AddValue("is_fixed",settings["is_fixed"])

        if settings.Has("fluid_pressure_type"):
            # if settings["hydrostatic"].GetBool() == False:
            if settings["fluid_pressure_type"].GetString() == "Uniform":
                self.params.AddValue("value",settings["value"])
                if settings["table"].GetInt() == 0:
                    self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.params)
                else:
                    self.params.AddValue("table",settings["table"])
                    self.process = KratosGeo.ApplyComponentTableProcess(self.model_part, self.params)
            elif settings["fluid_pressure_type"].GetString() == "Hydrostatic":
                self.params.AddValue("gravity_direction",settings["gravity_direction"])
                self.params.AddValue("reference_coordinate",settings["reference_coordinate"])
                self.params.AddValue("specific_weight",settings["specific_weight"])

                if settings.Has("pressure_tension_cut_off"):
                    self.params.AddValue("pressure_tension_cut_off",settings["pressure_tension_cut_off"])

                if settings.Has("is_seepage"):
                    self.params.AddValue("is_seepage",settings["is_seepage"])

                if settings["table"].GetInt() == 0:
                    self.process = KratosGeo.ApplyConstantHydrostaticPressureProcess(self.model_part, self.params)
                else:
                    self.params.AddValue("table",settings["table"])
                    self.process = KratosGeo.ApplyHydrostaticPressureTableProcess(self.model_part, self.params)
            elif settings["fluid_pressure_type"].GetString() == "Phreatic_Line":
                self.params.AddValue("gravity_direction",settings["gravity_direction"])
                self.params.AddValue("out_of_plane_direction",settings["out_of_plane_direction"])
                self.params.AddValue("first_reference_coordinate",settings["first_reference_coordinate"])
                self.params.AddValue("second_reference_coordinate",settings["second_reference_coordinate"])
                self.params.AddValue("specific_weight",settings["specific_weight"])

                if settings.Has("pressure_tension_cut_off"):
                    self.params.AddValue("pressure_tension_cut_off",settings["pressure_tension_cut_off"])

                if settings.Has("is_seepage"):
                    self.params.AddValue("is_seepage",settings["is_seepage"])

                if settings["table"][0].GetInt() == 0 and settings["table"][1].GetInt() == 0:
                    self.process = KratosGeo.ApplyConstantPhreaticLinePressureProcess(self.model_part, self.params)
                else:
                    self.params.AddValue("table",settings["table"])
                    self.process = KratosGeo.ApplyPhreaticLinePressureTableProcess(self.model_part, self.params)
            elif settings["fluid_pressure_type"].GetString() == "Interpolate_Line":
                self.params.AddValue("gravity_direction",settings["gravity_direction"])
                self.params.AddValue("out_of_plane_direction",settings["out_of_plane_direction"])

                if settings.Has("pressure_tension_cut_off"):
                    self.params.AddValue("pressure_tension_cut_off",settings["pressure_tension_cut_off"])

                if settings.Has("is_seepage"):
                    self.params.AddValue("is_seepage",settings["is_seepage"])

                if settings["table"].GetInt() == 0:
                    self.process = KratosGeo.ApplyConstantInterpolateLinePressureProcess(self.model_part, self.params)
                else:
                    self.process = KratosGeo.ApplyTimeDependentInterpolateLinePressureProcess(self.model_part, self.params)
            elif settings["fluid_pressure_type"].GetString() == "Phreatic_Surface":
                self.params.AddValue("gravity_direction",settings["gravity_direction"])
                self.params.AddValue("first_reference_coordinate",settings["first_reference_coordinate"])
                self.params.AddValue("second_reference_coordinate",settings["second_reference_coordinate"])
                self.params.AddValue("third_reference_coordinate",settings["third_reference_coordinate"])
                self.params.AddValue("specific_weight",settings["specific_weight"])

                if settings.Has("pressure_tension_cut_off"):
                    self.params.AddValue("pressure_tension_cut_off",settings["pressure_tension_cut_off"])

                if settings.Has("is_seepage"):
                    self.params.AddValue("is_seepage",settings["is_seepage"])

                if settings["table"][0].GetInt() == 0 and settings["table"][1].GetInt() == 0 and settings["table"][2].GetInt() == 0:
                    self.process = KratosGeo.ApplyConstantPhreaticSurfacePressureProcess(self.model_part, self.params)
                else:
                    self.params.AddValue("table",settings["table"])
                    self.process = KratosGeo.ApplyPhreaticSurfacePressureTableProcess(self.model_part, self.params)
            else:
                raise Exception("unkown fluid_pressure_type!")


        else:
            self.params.AddValue("value",settings["value"])
            if settings["table"].GetInt() == 0:
                self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.params)
            else:
                self.params.AddValue("table",settings["table"])
                self.process = KratosGeo.ApplyComponentTableProcess(self.model_part, self.params)

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
