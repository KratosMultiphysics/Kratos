# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.colors as colors

### This file contains functionalities that are commonly used in CoSimulation ###

def cs_print_info(label, *args):
    KM.Logger.PrintInfo(colors.bold(label), " ".join(map(str,args)))

def cs_print_warning(label, *args):
    KM.Logger.PrintWarning(colors.bold(label), " ".join(map(str,args)))


def SettingsTypeCheck(settings):
    if not isinstance(settings, KM.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
