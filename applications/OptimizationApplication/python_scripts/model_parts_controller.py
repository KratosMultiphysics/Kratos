# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KO

# ==============================================================================
def CreateController(model_parts_settings,model):
    return ModelPartsController(model_parts_settings,model)

# ==============================================================================
class ModelPartsController:
    # --------------------------------------------------------------------------
    def __init__(self, model_parts_settings,model):

        self.model_parts_settings = model_parts_settings
        self.model = model

        default_settings = KM.Parameters("""
        {
            "domain_size"           : 3,
            "model_part_name"       : "MODEL_PART_NAME",
            "model_import_settings"              : {
                "input_type"     : "mdpa",
                "input_filename" : "MODEL_PART_FILENAME"
            }
        }""")

      

        for itr in range(self.model_parts_settings.size()):
            for key in default_settings.keys():
                if not self.model_parts_settings[itr].Has(key):
                    raise RuntimeError("ModelPartsController: Required setting '{}' missing in 'model_part Nr.{}'!".format(key,itr+1))  
            self.model_parts_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["model_import_settings"].keys():
                if not self.model_parts_settings[itr]["model_import_settings"].Has(key):
                    raise RuntimeError("ModelPartsController: Required setting '{}' missing in 'model_import_settings' of 'model_part Nr.{}' !".format(key,itr+1))             
            self.model_parts_settings[itr]["model_import_settings"].ValidateAndAssignDefaults(default_settings["model_import_settings"]) 

        # here we initialize empty list of model parts, will be used in optmization
        self.model_parts = {}             

    # --------------------------------------------------------------------------
    def Initialize(self):
        self.__ImportModelParts()

    # --------------------------------------------------------------------------
    def CheckIfModelPartExists(self,model_part_name):
        exists = False
        for model_part_settings in self.model_parts_settings:
            if model_part_name == model_part_settings["model_part_name"].GetString():
                exists = True
                break
        return exists
        
    # --------------------------------------------------------------------------
    def SetMinimalBufferSize(self, buffer_size):
        for key,value in self.model_parts.items():
            if value.GetBufferSize() < buffer_size:
                value.SetBufferSize(buffer_size)
    # --------------------------------------------------------------------------
    def __ImportModelParts(self):

        for model_part_settings in self.model_parts_settings:
            model_part_name = model_part_settings["model_part_name"].GetString()
            input_type = model_part_settings["model_import_settings"]["input_type"].GetString()
            if input_type != "mdpa":
                raise RuntimeError("The model part "+format(model_part_name)+" for the optimization has to be read from the mdpa file!")

            optimization_model_part = self.model.CreateModelPart(model_part_name)
            optimization_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, model_part_settings["domain_size"].GetInt())

            input_filename = model_part_settings["model_import_settings"]["input_filename"].GetString()
            KM.ModelPartIO(input_filename).ReadModelPart(optimization_model_part)    
            self.model_parts[model_part_name] = optimization_model_part

        self.SetMinimalBufferSize(1)


# ==============================================================================
