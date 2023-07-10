import KratosMultiphysics

class DeleteModelPartsModeler(KratosMultiphysics.Modeler):

    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model = model
        self.settings = settings

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

        # Get the names of the model parts to be removed
        # Note that
        model_part_names = self.settings["model_part_names"].GetStringArray()
        model_part_names.sort()

        removed_sub_mp_list = []
        for model_part_name in model_part_names:
            if self.model.HasModelPart(model_part_name):
                # Save the name of the model part and sub model parts that are to be removed
                model_part = self.model.GetModelPart(model_part_name)
                aux_to_remove_sub_mp_list = self.__GetSubModelPartsNames(model_part)

                if self.settings["echo_level"].GetInt() > 1:
                    aux_info_msg = " and submodelparts:" if len(aux_to_remove_sub_mp_list) else "."
                    info_msg = f"About to delete model part '{model_part_name}'{aux_info_msg}"
                    for name in aux_to_remove_sub_mp_list:
                        info_msg += "\n - {}".format(name)
                    KratosMultiphysics.Logger.PrintInfo("DeleteModelPartModeler", info_msg)

                # Delete the model part through the model
                self.model.DeleteModelPart(model_part_name)

                if self.settings["echo_level"].GetInt() > 0:
                    KratosMultiphysics.Logger.PrintInfo("DeleteModelPartModeler", f"Model part '{model_part_name}' deleted.\n")

                # Accumulate the current deleted model parts to the
                removed_sub_mp_list += aux_to_remove_sub_mp_list

            elif model_part_name not in removed_sub_mp_list:
                warn_msg = "DeleteModelPartModeler", f"Asking to remove model part '{model_part_name}' but it is not present in model container."
                KratosMultiphysics.Logger.PrintWarning("DeleteModelPartModeler", warn_msg)

    def __GetDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters('''{
            "echo_level" : 0,
            "model_part_names" : [""]
        }''')
        return default_settings

    def __GetSubModelPartsNames(self, model_part):
        aux_list = []
        for sub_mp in model_part.SubModelParts:
            aux_list.append(sub_mp.FullName())
            aux_list += [name for name in self.__GetSubModelPartsNames(sub_mp)]
        return aux_list

def Factory(model, settings):
    return DeleteModelPartsModeler(model, settings)
