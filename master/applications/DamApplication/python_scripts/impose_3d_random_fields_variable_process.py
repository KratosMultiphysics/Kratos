import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

try:
    from gstools import SRF, Gaussian
except ImportError:
    raise ImportError("The use of the random fields module requires 'gstools'!")

from statistics import mean, variance
from math import sqrt

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return Impose3dRandomFieldsVariableProcess(Model, settings["Parameters"])

class Impose3dRandomFieldsVariableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        mean_value = settings["mean_value"].GetDouble()
        min_value = settings["min_value"].GetDouble()
        max_value = settings["max_value"].GetDouble()
        var = settings["variance"].GetDouble()
        corr_length = settings["corr_length"].GetInt()

        # Gaussian random field with exponential covariance
        model = Gaussian(dim = 3, var = var, len_scale = corr_length)
        srf = SRF(model)

        x = []
        y = []
        z = []
        ids = []

        for node in model_part.Nodes:
            x.append(node.X)
            y.append(node.Y)
            z.append(node.Z)
            ids.append(node.Id)

        field = srf((x, y, z))

        field_mean = mean(field)
        field_var = variance(field, field_mean)

        variable_values = []

        for field_i in field:
            variable_values.append(mean_value + ((field_i-field_mean)*var/sqrt(field_var)))

        ### Truncate values
        for variable_value in variable_values:
            if variable_value < min_value:
                variable_value = min_value
            if variable_value > max_value:
                variable_value = max_value

        self.table = KratosMultiphysics.PiecewiseLinearTable()

        for i in range(len(variable_values)):
            self.table.AddRow(int(ids[i]), float(variable_values[i]))

        self.process = KratosDam.DamRandomFieldsVariableProcess(model_part, self.table, settings)


    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
