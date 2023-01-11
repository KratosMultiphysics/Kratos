# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# importing the Kratos Library
import KratosMultiphysics as KM

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

        self.root_model_part_name_filename_map={}
        self.name_root_model_part_map={}
        for model_part_settings in self.model_parts_settings:
            model_part_name = model_part_settings["name"].GetString()
            input_filename = model_part_settings["settings"]["input_filename"].GetString()
            input_type = model_part_settings["type"].GetString()
            if input_type != "mdpa":
                raise RuntimeError("The model part "+format(model_part_name)+" for the optimization has to be read from the mdpa file!")
            if model_part_name in self.root_model_part_name_filename_map.keys():
                raise RuntimeError("ModelPartsController: there are duplicated model parts !")
            if input_filename in self.root_model_part_name_filename_map.values():
                raise RuntimeError("ModelPartsController: there are duplicated input files !")
            optimization_model_part = self.model.CreateModelPart(model_part_name)
            optimization_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, model_part_settings["settings"]["domain_size"].GetInt())
            self.root_model_part_name_filename_map[model_part_name] = input_filename
            self.name_root_model_part_map[model_part_name] = optimization_model_part

    # --------------------------------------------------------------------------
    def Initialize(self):
        self.__ImportRootModelParts()

    # --------------------------------------------------------------------------
    def CheckIfRootModelPartExists(self,root_model_part_name,raise_error=True):
        extracted_root_model_part_name = root_model_part_name.split(".")[0]
        if extracted_root_model_part_name in self.name_root_model_part_map.keys():
            return True
        else:
            if raise_error:
                raise RuntimeError("ModelPartsController: CheckIfRootModelPartExists: Root model part {} does not exist!".format(extracted_root_model_part_name))
            else:
                return False
    # --------------------------------------------------------------------------
    def CheckIfRootModelPartsExist(self,root_model_parts_name,raise_error=True):
        if type(root_model_parts_name) is not list:
            raise RuntimeError("ModelPartsController: CheckIfRootModelPartsExist requires list of model parts")

        if_exist = True
        for root_model_part_name in root_model_parts_name:
            extracted_root_model_part_name = root_model_part_name.split(".")[0]
            if not extracted_root_model_part_name in self.name_root_model_part_map.keys():
                if raise_error:
                    raise RuntimeError("ModelPartsController: CheckIfRootModelPartsExist: Root model part {} does not exist!".format(extracted_root_model_part_name))
                else:
                    if_exist = False
                    break

        return if_exist
    # --------------------------------------------------------------------------
    def GetRootModelPart(self, root_model_part_name):
        extracted_root_model_part_name = root_model_part_name.split(".")[0]
        if not extracted_root_model_part_name in self.name_root_model_part_map.keys():
            raise RuntimeError("ModelPartsController: Try to get root model part {} which does not exist.".format(root_model_part_name))
        else:
            return self.name_root_model_part_map[extracted_root_model_part_name]
    # --------------------------------------------------------------------------
    def UpdateTimeStep(self, step):
        for root_model_part in self.name_root_model_part_map.values():
            root_model_part.CloneTimeStep(step)
            root_model_part.ProcessInfo.SetValue(KM.STEP, step)
    # --------------------------------------------------------------------------
    def SetMinimalBufferSize(self, buffer_size):
        for root_model in self.name_root_model_part_map.values():
            if root_model.GetBufferSize() < buffer_size:
                root_model.SetBufferSize(buffer_size)
    # --------------------------------------------------------------------------
    def __ImportRootModelParts(self):

        for name,root_model in self.name_root_model_part_map.items():
            input_filename = self.root_model_part_name_filename_map[name]
            KM.ModelPartIO(input_filename).ReadModelPart(root_model)

        self.SetMinimalBufferSize(1)


# ==============================================================================
