import KratosMultiphysics as KM
import kerastuner as kt
import numpy as np


class LoadTunableParametersUtility:

    def __init__(self, settings):

        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing layer settings.
        """

        default_settings = KM.Parameters("""{
            "type_of_tuning"           : "",
            "name"                     : "",
            "values"                   : [],
            "min_value"                : 0.0,
            "max_value"                : 0.0,
            "step"                     : 0.0,
            "sampling"                 : ""
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.type_of_tuning = settings["type_of_tuning"].GetString()
        self.name = settings["name"].GetString()

        if self.type_of_tuning == "Choice":
            if settings["values"].IsVector():
                self.values = []
                values_vector = settings["values"].GetVector()
                for i in range(values_vector.Size()):
                    self.values.append(values_vector[i])
            else:
                self.values = settings["values"].GetStringArray()

        if self.type_of_tuning == "Int":
            self.min_value = settings["min_value"].GetInt()
            self.max_value = settings["max_value"].GetInt()
            self.step = settings["step"].GetInt()
            if not self.step > 1:
                self.step = 1
            self.sampling = settings["sampling"].GetString()
            if self.sampling == "":
                self.sampling = None

        if self.type_of_tuning == "Float":
            self.min_value = settings["min_value"].GetDouble()
            self.max_value = settings["max_value"].GetDouble()
            self.step = settings["step"].GetDouble()
            if not self.step > 0.0:
                self.step = None
            self.sampling = settings["sampling"].GetString()
            if self.sampling == "":
                self.sampling = None

    def Load(self, hp):
        """ Loads the settings for the tuning of hyperparameter hp """
        
        if self.type_of_tuning == "Boolean":
            return hp.Boolean(self.name)

        elif self.type_of_tuning == "Choice":
            return hp.Choice(self.name, values = self.values)

        elif self.type_of_tuning == "Int":
            return hp.Int(self.name, self.min_value, self.max_value, step = self.step, sampling = self.sampling)

        elif self.type_of_tuning == "Float":
            return hp.Float(self.name, self.min_value, self.max_value, step = self.step, sampling = self.sampling)
        else:
            raise Exception(self.type_of_tuning + " is not a possible type of tuning. Available types are Boolean, Choice, Int and Float.")