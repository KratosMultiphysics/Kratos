import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, Logger

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.PlasmaDynamicsApplication as PlasmaDynamics
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from KratosMultiphysics.SwimmingDEMApplication.variables_management import VariablesManager as SDEMVariablesManager


def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)


def GetGlobalVariableByName(variable_name):
    modules = [Kratos, DEM, PlasmaDynamics]
    for mod in modules:
        try:
            return Kratos.KratosGlobals.GetVariable(variable_name)
        except Exception:
            pass
    names = [mod.__name__ for mod in modules]
    error_message = ('No variable with name \''
        + variable_name + '\' exists in either of the modules:\n')
    for name in names[:-1]:
        error_message += name + ', '
    error_message += 'or ' + names[-1] + '.'
    raise AttributeError(error_message)



#TODO: complete the class
class VariablesManager(SDEMVariablesManager):
    
    def __init__(self, parameters):
        super.__init__(parameters)

    
