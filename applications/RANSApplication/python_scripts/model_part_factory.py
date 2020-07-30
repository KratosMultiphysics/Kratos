import KratosMultiphysics as Kratos
from KratosMultiphysics.RANSApplication import RansVariableUtilities


def CreateDuplicateModelPart(origin_modelpart, destination_modelpart_name,
                             element_name, condition_name,
                             original_condition_name):
    # domain_size = origin_modelpart.ProcessInfo[Kratos.DOMAIN_SIZE]
    model = origin_modelpart.GetModel()
    connectivity_preserve_modeler = Kratos.ConnectivityPreserveModeler()

    if not model.HasModelPart(destination_modelpart_name):
        model_part = model.CreateModelPart(destination_modelpart_name)

        # TODO: Remove this line once the warnings from connectivity preserve modeller is gone, otherwise,
        #       output will be cluttered with lots of missing variable warnings
        RansVariableUtilities.CopyNodalSolutionStepVariablesList(
            origin_modelpart, model_part)

        # TODO: [PeriodicCondition]
        #       Currently, all the conditions will be replaced with the given new condition. This is an issue
        #       in the case of periodic cases in mpi, there we have to put PeriodicConditions in the mdpa file,
        #       where MetisParitioner will use that condition list to properly partition it. Therefore, "original_condition_name"
        #       is not used in this method at the moment.
        #       Following is one of the proposals to make PeriodicConditions to work with connectivity_preserve_modeller.
        # connectivity_preserve_modeler.GenerateModelPart(
        #     origin_modelpart, model_part, element_name, condition_name,
        #     original_condition_name + str(domain_size) + "D" + str(domain_size)
        #     + "N")

        connectivity_preserve_modeler.GenerateModelPart(
            origin_modelpart, model_part, element_name, condition_name)

    Kratos.Logger.PrintInfo("RANSModelPartFactory",
                            "Created " + destination_modelpart_name)
    return model.GetModelPart(destination_modelpart_name)


def ApplyFlagsToNodeList(nodes, flag_name, value):
    allowed_flag_names = ["inlet", "outlet", "structure", "none"]

    if not flag_name in allowed_flag_names:
        raise Exception("Unknown flag name: " + flag_name +
                        ". Allowed flag names are: " + allowed_flag_names)

    variable_utils = Kratos.VariableUtils()
    if (flag_name == "inlet"):
        variable_utils.SetFlag(Kratos.INLET, value, nodes)
    elif (flag_name == "outlet"):
        variable_utils.SetFlag(Kratos.OUTLET, value, nodes)
    elif (flag_name == "structure"):
        variable_utils.SetFlag(Kratos.STRUCTURE, value, nodes)


def ApplyFlagsToConditionsList(conditions, flag_name, value):
    allowed_flag_names = ["inlet", "outlet", "structure", "none"]

    if not flag_name in allowed_flag_names:
        raise Exception("Unknown flag name: " + flag_name +
                        ". Allowed flag names are: " + allowed_flag_names)

    variable_utils = Kratos.VariableUtils()
    if (flag_name == "inlet"):
        variable_utils.SetFlag(Kratos.INLET, value, conditions)
    elif (flag_name == "outlet"):
        variable_utils.SetFlag(Kratos.OUTLET, value, conditions)
    elif (flag_name == "structure"):
        variable_utils.SetFlag(Kratos.STRUCTURE, value, conditions)
