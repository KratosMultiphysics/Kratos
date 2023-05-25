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
from KratosMultiphysics.OptimizationApplication.model_part_controllers.mdpa_model_part_controller import MdpaModelPartController

# ==============================================================================
def CreateController(model_parts_settings: KM.Parameters, model: KM.Model):
    return ModelPartsController(model_parts_settings, model)

# ==============================================================================
class ModelPartsController:
    # --------------------------------------------------------------------------
    def __init__(self, model_parts_settings: KM.Parameters, model: KM.Model):

        self.model_parts_settings = model_parts_settings
        self.model = model

        self.mdpa_model_part_controllers: 'list[MdpaModelPartController]' = []
        for params in self.model_parts_settings:
            params["settings"].AddString("model_part_name", params["name"].GetString())
            self.mdpa_model_part_controllers.append(MdpaModelPartController(model, params["settings"]))

    # --------------------------------------------------------------------------
    def Initialize(self) -> None:
        for mdpa_model_part_controller in self.mdpa_model_part_controllers:
            mdpa_model_part_controller.ImportModelPart()

    # --------------------------------------------------------------------------
    def CheckIfRootModelPartsExist(self, root_model_parts_name: 'list[str]', raise_error = True) -> bool:
        if not isinstance(root_model_parts_name, list):
            raise RuntimeError("ModelPartsController: CheckIfRootModelPartsExist requires list of model parts")

        if_exist = True
        for root_model_part_name in root_model_parts_name:
            extracted_root_model_part_name = root_model_part_name.split(".")[0]
            if not self.model.HasModelPart(extracted_root_model_part_name):
                if raise_error:
                    raise RuntimeError("ModelPartsController: CheckIfRootModelPartsExist: Root model part {} does not exist!".format(extracted_root_model_part_name))
                else:
                    if_exist = False
                    break

        return if_exist
    # --------------------------------------------------------------------------
    def GetModelPart(self, model_part_name: str) -> KM.ModelPart:
        if not model_part_name in self.model.GetModelPartNames():
            raise RuntimeError("ModelPartsController: Try to get model part {} which does not exist.".format(model_part_name))
        else:
            return self.model.GetModelPart(model_part_name)
    # --------------------------------------------------------------------------
    def GetRootModelPart(self, root_model_part_name: str) -> KM.ModelPart:
        extracted_root_model_part_name = root_model_part_name.split(".")[0]
        if not self.model.HasModelPart(extracted_root_model_part_name):
            raise RuntimeError("ModelPartsController: Try to get root model part {} which does not exist.".format(root_model_part_name))
        else:
            return self.model[extracted_root_model_part_name]
    # --------------------------------------------------------------------------
    def GetRootModelParts(self, root_model_parts_name: 'list[str]') -> 'list[KM.ModelPart]':
        if not isinstance(root_model_parts_name, list):
            raise RuntimeError("ModelPartsController: GetRootModelParts requires list of model parts name")

        list_root_model_parts = []
        for root_model_part_name in root_model_parts_name:
            extracted_root_model_part_name = root_model_part_name.split(".")[0]
            if not self.model.HasModelPart(extracted_root_model_part_name):
                raise RuntimeError("ModelPartsController: GetRootModelParts: Root model part {} does not exist!".format(extracted_root_model_part_name))
            else:
                list_root_model_parts.append(self.model[extracted_root_model_part_name])

        return list_root_model_parts

    # --------------------------------------------------------------------------
    def UpdateTimeStep(self, step: int) -> None:
        for mdpa_odel_part_controller in self.mdpa_model_part_controllers:
            mdpa_odel_part_controller.GetModelPart().CloneTimeStep(step)
            mdpa_odel_part_controller.GetModelPart().ProcessInfo.SetValue(KM.STEP, step)
    # --------------------------------------------------------------------------
    def SetMinimalBufferSize(self, buffer_size: int) -> None:
        for mdpa_odel_part_controller in self.mdpa_model_part_controllers:
            if mdpa_odel_part_controller.GetModelPart().GetBufferSize() < buffer_size:
                mdpa_odel_part_controller.GetModelPart().SetBufferSize(buffer_size)

# ==============================================================================
