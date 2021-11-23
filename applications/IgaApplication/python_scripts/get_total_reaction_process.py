import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SM

from matplotlib import pyplot as plt
import matplotlib.tri as tri
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.animation as animation
import numpy as np

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OuputTotalReactionProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class OuputTotalReactionProcess(KM.Process):
    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "model_part_name" : "please_specify_model_part_name",
            "output_file_name" : "total_reactions",
            "condition_variable_list" : []
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

        self.output_file_name = settings["output_file_name"].GetString() + ".txt"
        with open(self.output_file_name, 'w') as output_file:
            output_file.write("")

        self.condition_results_scalar, self.condition_results_vector = \
            CreateVariablesListFromInput(settings["condition_variable_list"])

        self.time_data = []

    def ExecuteFinalizeSolutionStep(self):
        time_step = self.model_part.ProcessInfo[KM.TIME]
        self.time_data.append(time_step)

        with open(self.output_file_name, 'a') as output_file:
            output_file.write(str(time_step))

            for variable in self.condition_results_scalar:
                value = 0.0
                for condition in self.model_part.Conditions:
                    value += condition.CalculateOnIntegrationPoints(variable, self.model_part.ProcessInfo)[0]
                output_file.write("  " + str(value))

            for variable in self.condition_results_vector:
                value_x = 0.0
                value_y = 0.0
                value_z = 0.0
                value = 0.0
                for condition in self.model_part.Conditions:
                    result = condition.CalculateOnIntegrationPoints(variable, self.model_part.ProcessInfo)[0]
                    value_x += result[0]
                    value_y += result[1]
                    value_z += result[2]
                    value += sqrt(result[0]**2 + result[1]**2 + result[2]**2)
                output_file.write("  " + str(value_x) + "  " + str(value_y) + "  " + str(value_z) + "  " + str(value))

            output_file.write("\n")

def CreateVariablesListFromInput(param):
    '''Parse a list of variables from input.'''
    scalar_variables = []
    vector_variables = []
    admissible_scalar_types = ["Bool", "Integer", "Unsigned Integer", "Double"]
    admissible_vector_types = ["Array", "Vector"]

    variable_list = KM.kratos_utilities.GenerateVariableListFromInput(param)

    # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
    for variable in variable_list:
        if KM.KratosGlobals.GetVariableType(variable.Name()) in admissible_scalar_types:
            scalar_variables.append(variable)
        elif KM.KratosGlobals.GetVariableType(variable.Name()) in admissible_vector_types:
            vector_variables.append(variable)
        else:
            raise Exception("unsupported variables type: " + str(type(variable)))

    return scalar_variables, vector_variables