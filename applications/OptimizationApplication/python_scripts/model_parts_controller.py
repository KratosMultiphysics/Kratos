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
            "name"       : "MODEL_PART_NAME",
            "type"     : "mdpa",
            "settings"              : {
                "domain_size"           : 3,
                "input_filename" : "MODEL_PART_FILENAME"
            }
        }""")      

        for itr in range(self.model_parts_settings.size()):
            for key in default_settings.keys():
                if not self.model_parts_settings[itr].Has(key):
                    raise RuntimeError("ModelPartsController: Required setting '{}' missing in 'model_part Nr.{}'!".format(key,itr+1))  
            self.model_parts_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["settings"].keys():
                if not self.model_parts_settings[itr]["settings"].Has(key):
                    raise RuntimeError("ModelPartsController: Required setting '{}' missing in 'model_import_settings' of 'model_part Nr.{}' !".format(key,itr+1))             
            self.model_parts_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"]) 
  
        self.model_parts_name=[]
        self.model_parts_input_filename=[]
        for model_part_settings in self.model_parts_settings:
            model_part_name = model_part_settings["name"].GetString()
            input_filename = model_part_settings["settings"]["input_filename"].GetString()
            input_type = model_part_settings["type"].GetString()
            if input_type != "mdpa":
                raise RuntimeError("The model part "+format(model_part_name)+" for the optimization has to be read from the mdpa file!")
            if model_part_name in self.model_parts_name:
                raise RuntimeError("ModelPartsController: there are duplicated model parts !") 
            if input_filename in self.model_parts_input_filename:
                raise RuntimeError("ModelPartsController: there are duplicated input files !")                
            optimization_model_part = self.model.CreateModelPart(model_part_name)
            optimization_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, model_part_settings["settings"]["domain_size"].GetInt())   
            self.model_parts_name.append(model_part_name) 
            self.model_parts_input_filename.append(input_filename)              

    # --------------------------------------------------------------------------
    def Initialize(self):
        self.__ImportModelParts()

    # --------------------------------------------------------------------------
    def CheckIfModelPartExists(self,model_part_name):
        exists = False
        for model_part_settings in self.model_parts_settings:
            if model_part_name == model_part_settings["name"].GetString():
                exists = True
                break
        return exists

    # --------------------------------------------------------------------------
    def GetModelPart(self, model_part_name):
        if not model_part_name in self.model.GetModelPartNames():
            raise RuntimeError("AnalyzersController: Try to get model part {} which does not exist.".format(model_part_name))
        else:
            return self.model.GetModelPart(model_part_name)

    # --------------------------------------------------------------------------
    def SetMinimalBufferSize(self, buffer_size):
        for model_part_name in self.model.GetModelPartNames():
            optimization_model_part = self.model.GetModelPart(model_part_name)
            if optimization_model_part.GetBufferSize() < buffer_size:
                optimization_model_part.SetBufferSize(buffer_size)
    # --------------------------------------------------------------------------
    def __ImportModelParts(self):

        for model_part_settings in self.model_parts_settings:
            model_part_name = model_part_settings["name"].GetString()
            optimization_model_part = self.model.GetModelPart(model_part_name)
            input_filename = model_part_settings["settings"]["input_filename"].GetString()
            KM.ModelPartIO(input_filename).ReadModelPart(optimization_model_part)

        self.SetMinimalBufferSize(1)


# ==============================================================================
