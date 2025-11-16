import KratosMultiphysics
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning


def Factory(settings, model):
    if (not isinstance(model, KratosMultiphysics.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return ApplyWallFunctionProcess(model, settings["Parameters"])


class ApplyWallFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        """This process sets wall function parameters

        This process is used to activate wall function in wall conditions by setting
        RANS_IS_WALL_FUNCTION_ACTIVE to 1 in conditions. Additionally it also sets STRUCTURE and SLIP flag
        on the same nodes as well as conditions in same model part.

        Previously, all the wall conditions should be in one model part/sub model part. Having wall conditions
        in more than one model part/sub model part is not allowed. This has been deprecated.

        Currently, the wall conditions can be in different modelpart/submodelpart and should be passed as a
        list.

        Args:
            Model (Kratos.Model): [description]
            settings (Kratos.Parameters): [description]

        Raises:
            Exception: [description]
        """

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_names" : [],
                "all_wall_model_part_name": "ALL_WALL_MODEL_PART",
                "activate_wall_functions" : true
            }  """)

        #  backward compatibility: "model_part_name" --> "model_part_names"
        if settings.Has("model_part_name"):
            if settings.Has("model_part_names"):
                if settings["model_part_names"].size() > 0:
                    err_msg = "The parameters in {} cannot have both \n".format(self.__class__.__name__)
                    err_msg += "'model_part_name' and 'model_part_names'.\n'model_part_name' is deprecated "
                    err_msg += "and is only allowed for backward compatibility. "
                    err_msg += "It is recommended to use just 'model_part_names'.\n"
                    err_msg += "{}".format(settings)
                    raise Exception(err_msg)

            else:
                # if settings has only "model_part_name", issue a deprecation warning.
                warn_msg = "\n'model_part_name' is deprecated and is only allowed for backward compatibility."
                warn_msg += "\nIt is recommended to use 'model_part_names' instead."
                warn_msg += "\n{}".format(settings)
                # IssueDeprecationWarning(str(self.__class__.__name__), warn_msg)

                # add "model_part_names" and remove "model_part_name" from settings
                settings.AddEmptyList("model_part_names")
                settings["model_part_names"].Append(settings["model_part_name"])
                settings.RemoveValue("model_part_name")

        settings.ValidateAndAssignDefaults(default_parameters)

        if settings["model_part_names"].size() == 0:
            err_msg = "Empty 'model_part_names' parameter in {}.\n".format(self.__class__.__name__)
            err_msg += "Must contain atleast one wall modelpart.\n"
            err_msg += "{}".format(settings)
            raise Exception(err_msg)

        self.wall_mp_list = settings["model_part_names"].GetStringArray()
        self.fluid_model_part = model[self.wall_mp_list[0]].GetRootModelPart()

        if (not self.fluid_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]):
            # TODO: Remove this block. This is only required for time averaged adjoint analysis primal runs only with turbulence models
            # this is a hack for the time being
            if not self.fluid_model_part.ProcessInfo.Has(KratosRANS.WALL_MODEL_PART_NAME):
                all_wall_model_part_name = settings["all_wall_model_part_name"].GetString()
                if all_wall_model_part_name.find(".") == -1:
                    wall_mp = self.fluid_model_part.CreateSubModelPart(settings["all_wall_model_part_name"].GetString())
                else:
                    wall_mp = model.CreateModelPart(all_wall_model_part_name)
                self.fluid_model_part.ProcessInfo.SetValue(KratosRANS.WALL_MODEL_PART_NAME, wall_mp.FullName())

            self.model_part = model[self.fluid_model_part.ProcessInfo[KratosRANS.WALL_MODEL_PART_NAME]]
            if (self.model_part.NumberOfNodes() != 0):
                raise Exception(
                    "ApplyWallFunctionProcess can be applied only once. Therefore please group all wall model parts to one main model part and apply this process to it."
                )

            # Create a submodel part that stores all the wall model parts provided in the list
            self.wall_smp_name = self.model_part.FullName()
            self.fluid_model_part.ProcessInfo.SetValue(KratosRANS.WALL_MODEL_PART_NAME, self.wall_smp_name)

            # Add each wall from the list to the ALL_WALL submodelpart
            for wall_mp_name in self.wall_mp_list:
                # Get the wall model part in the list
                wall_modelpart = model[wall_mp_name]

                # making list containing node, element and conditions IDs respectively of each wallmodelpart
                # [inspired by <chimera_modelpart_import.py> script]
                wmp_node_id_array       = [node.Id for node in wall_modelpart.Nodes]
                wmp_element_id_array    = [element.Id for element in wall_modelpart.Elements]
                wmp_condition_id_array  = [condition.Id for condition in wall_modelpart.Conditions]

                # Add the entities of wallmodelpart into the ALL_WALL submodelpart
                self.model_part.AddNodes(wmp_node_id_array)
                self.model_part.AddElements(wmp_element_id_array)
                self.model_part.AddConditions(wmp_condition_id_array)

            for node in self.model_part.Nodes:
                node.Set(KratosMultiphysics.STRUCTURE, True)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.STRUCTURE, True)

            self.model_part.SetValue(KratosRANS.RANS_IS_WALL_FUNCTION_ACTIVE, settings["activate_wall_functions"].GetBool())

            if (settings["activate_wall_functions"].GetBool()):
                for node in self.model_part.Nodes:
                    node.Set(KratosMultiphysics.SLIP, True)
                    node.SetValue(KratosMultiphysics.Y_WALL, 1.0)

                for condition in self.model_part.Conditions:
                    condition.Set(KratosMultiphysics.SLIP, True)
                    condition.SetValue(KratosRANS.RANS_IS_WALL_FUNCTION_ACTIVE, 1)
            else:
                for condition in self.model_part.Conditions:
                    condition.SetValue(KratosRANS.RANS_IS_WALL_FUNCTION_ACTIVE, 0)
