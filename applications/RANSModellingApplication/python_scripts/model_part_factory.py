import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSModellingApplication as KratosRANS


def CreateDuplicateModelPart(origin_modelpart, destination_modelpart_name,
                             element_name, condition_name):
    domain_size = origin_modelpart.ProcessInfo[Kratos.DOMAIN_SIZE]
    model = origin_modelpart.GetModel()
    connectivity_preserve_modeler = Kratos.ConnectivityPreserveModeler()

    if not model.HasModelPart(destination_modelpart_name):
        model_part = model.CreateModelPart(destination_modelpart_name)
        KratosRANS.RansVariableUtils().CopyNodalSolutionStepVariablesList(
            origin_modelpart, model_part)

        connectivity_preserve_modeler.GenerateModelPart(
            origin_modelpart, model_part, element_name, condition_name)

    return model.GetModelPart(destination_modelpart_name)


def ApplyFlagsToModelPart(model_part, flag_name, value):
    allowed_flag_names = ["inlet", "outlet", "structure"]

    if not flag_name in allowed_flag_names:
        raise Exception("Unknown flag name: " + flag_name + ". Allowed flag names are: " + allowed_flag_names)

    variable_utils = Kratos.VariableUtils()
    if (flag_name == "inlet"):
        variable_utils.SetFlag(Kratos.INLET, value, model_part.Nodes)
    elif (flag_name == "outlet"):
        variable_utils.SetFlag(Kratos.OUTLET, value, model_part.Nodes)
    elif (flag_name == "structure"):
        variable_utils.SetFlag(Kratos.STRUCTURE, value, model_part.Nodes)
