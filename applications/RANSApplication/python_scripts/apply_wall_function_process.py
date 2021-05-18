import KratosMultiphysics
import KratosMultiphysics.RANSApplication as KratosRANS


def Factory(settings, Model):
    if (not isinstance(Model, KratosMultiphysics.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return ApplyWallFunctionProcess(Model, settings["Parameters"])


class ApplyWallFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        """This process sets wall function parameters

        This process is used to activate wall function in wall conditions by setting
        RANS_IS_WALL_FUNCTION_ACTIVE to 1 in conditions. Additionally it also sets STRUCTURE and SLIP flag
        on the same nodes as well as conditions in same model part.

        Currently, all the wall conditions should be in one model part/sub model part. Having wall conditions
        in more than one model part/sub model part is not allowed.

        Args:
            Model (Kratos.Model): [description]
            settings (Kratos.Parameters): [description]

        Raises:
            Exception: [description]
        """

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME"
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        if (not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]):            
            process_info = self.model_part.ProcessInfo
            if (process_info.Has(KratosRANS.WALL_MODEL_PART_NAME)):
                raise Exception(
                    "ApplyWallFunctionProcess can be applied only once. Therefore please group all wall model parts to one main model part and apply this process to it."
                )

            self.model_part.ProcessInfo.SetValue(
                KratosRANS.WALL_MODEL_PART_NAME,
                settings["model_part_name"].GetString())

            for node in self.model_part.Nodes:
                node.Set(KratosMultiphysics.SLIP, True)
                node.Set(KratosMultiphysics.STRUCTURE, True)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.SLIP, True)
                condition.Set(KratosMultiphysics.STRUCTURE, True)
                condition.SetValue(KratosRANS.RANS_IS_WALL_FUNCTION_ACTIVE, 1)
