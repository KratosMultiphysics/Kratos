from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from importlib import import_module

################# Please do not change the following class:
class KratosModelerFactory(object):
    def ConstructListOfModelers( self, modeler_list ):
        constructed_modelers = []
        for modeler_item in modeler_list:
            if not modeler_item.Has("python_module"):
                KM.Logger.PrintWarning("Your list of modelers: ", modeler_list)
                raise NameError('"python_module" must be defined in your parameters. Check all your modelers')

            # python-script that contains the modeler
            python_module_name = modeler_item["python_module"].GetString()

            if modeler_item.Has("kratos_module"): # for Kratos-Modelers
                """Location of the modelers in Kratos; eg.:
                - KratosMultiphysics
                - KratosMultiphysics.FluidDynamicsApplication
                """

                kratos_module_name = modeler_item["kratos_module"].GetString()
                if not kratos_module_name.startswith("KratosMultiphysics"):
                    kratos_module_name = "KratosMultiphysics." + kratos_module_name

                full_module_name = kratos_module_name + "." + python_module_name
                python_module = import_module(full_module_name)

                p = python_module.Factory(modeler_item)
                constructed_modelers.append( p )

            else: # for user-defined modelers
                python_module = import_module(python_module_name)
                p = python_module.Factory(modeler_item)
                constructed_modelers.append( p )

        return constructed_modelers


################# C++ MODELERS ###################
def Factory(settings):
    return ModelerFactory.Create(settings["modeler_name"].GetString(), settings["Parameters"])
