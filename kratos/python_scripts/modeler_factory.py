# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from importlib import import_module

class KratosModelerFactory(object):
    def ConstructListOfModelers( self, model, modeler_list ):
        constructed_modelers = []
        for modeler_item in modeler_list:
            if modeler_item.Has("modeler_name"):
                modeler_name = modeler_item["modeler_name"].GetString()
                if (KM.HasModeler(modeler_name)):
                    constructed_modelers.append( KM.CreateModeler(modeler_name, model, modeler_item["Parameters"]) )

                else:
                    if modeler_item.Has("kratos_module"):
                        """Location of the modelers in Kratos; eg.:
                        - KratosMultiphysics
                        - KratosMultiphysics.FluidDynamicsApplication
                        """
                        kratos_module_name = modeler_item["kratos_module"].GetString()

                        if not kratos_module_name.startswith("KratosMultiphysics"):
                            kratos_module_name = "KratosMultiphysics." + kratos_module_name
                        modeler_name = kratos_module_name + ".modelers." + modeler_name


                    python_module = import_module(modeler_name)
                    modeler = python_module.Factory(model, modeler_item["Parameters"])
                    constructed_modelers.append(modeler)

        return constructed_modelers
