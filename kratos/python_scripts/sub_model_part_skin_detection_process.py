# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    main_model_part = Model.GetModelPart(
        settings["Parameters"]["selection_settings"]["sub_model_part_name"].GetString()
        ).GetRootModelPart()
    domain_size = main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

    if domain_size == 2:
        return KratosMultiphysics.SubModelPartSkinDetectionProcess2D(main_model_part,settings["Parameters"])
    else:
        return KratosMultiphysics.SubModelPartSkinDetectionProcess3D(main_model_part,settings["Parameters"])
