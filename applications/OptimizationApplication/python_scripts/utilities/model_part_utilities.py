from typing import Any
from enum import Enum
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class ModelPartOperation:
    class OperationType(Enum):
        UNION = 1,
        INTERSECT = 2

    def __init__(self, model: Kratos.Model, operation_type: OperationType, suggested_model_part_name: str, list_of_operation_model_part_full_names: 'list[str]', add_neighbours: bool) -> None:
        self.model = model
        self.operation_type = operation_type
        self.suggested_model_part_name = suggested_model_part_name
        self.list_of_operation_model_part_full_names = sorted(list(set(list_of_operation_model_part_full_names)))
        self.add_neighbours = add_neighbours

        if len(list_of_operation_model_part_full_names) == 0:
            raise RuntimeError("No operating model part names are provided.")

        if suggested_model_part_name.find("#") != -1 or any([mp_name.find("#") != -1 for mp_name in list_of_operation_model_part_full_names]):
            # find '#' in the model part names which is not allwed.
            raise RuntimeError(f"The model part names cannot contain '#' character. Parsed model part names are as followings:\n\t suggested model part name: \"{suggested_model_part_name}\"\n\toperation model part names:\n\t" + "\n\t\t".join(list_of_operation_model_part_full_names))

        # set the root model part. This needs to always exist.
        self.root_model_part = self.model[self.list_of_operation_model_part_full_names[0].split(".")[0]]
        self.status_msg = ""

    def GetRootModelPart(self) -> Kratos.ModelPart:
        return self.root_model_part

    def GetModelPartFullName(self) -> str:
        # check whether is there a need to perform binary operations
        if self.operation_type == ModelPartOperation.OperationType.UNION or self.operation_type == ModelPartOperation.OperationType.INTERSECT:
            if len(self.list_of_operation_model_part_full_names) == 1:
                # all the operation model parts are the same, hence no need of doing an operation.
                self.status_msg = ""
                return self.list_of_operation_model_part_full_names[0]

        # now check in the status messages of the root model part whether this operation is already added.
        status_msg_prefix = "ModelPartUtilities_created#"
        status_msg_suffix = f"#{self.operation_type.name}#{self.add_neighbours:d}#" + "#".join(self.list_of_operation_model_part_full_names)
        status_msg_log = KratosOA.OptAppModelPartUtils.GetModelPartStatusLog(self.GetRootModelPart())
        for status_msg in status_msg_log:
            if status_msg.startswith(status_msg_prefix) and status_msg.endswith(status_msg_suffix):
                # found the same operation done on a different model part, hence sending that name
                self.status_msg = status_msg
                full_model_part_name = f"{self.GetRootModelPart().FullName()}.{status_msg.split('#')[1]}"
                Kratos.Logger.PrintInfo("ModelPartUtilities", f"Using \"{full_model_part_name}\" with the same operation instead of suggested model part with name = \"{self.GetRootModelPart().FullName()}.{self.suggested_model_part_name}\".")
                return full_model_part_name

        # it is not already called for creation. Then put that in the status msg.
        if self.root_model_part.HasSubModelPart(self.suggested_model_part_name):
            # this means, it already has a model part with suggeted name, but
            # it does not match the operation identifier. So throw an error
            raise RuntimeError(f"Found an already existing submodel part named \"{self.suggested_model_part_name}\" in {self.root_model_part.FullName()} without the required operation identifier = \"{status_msg_suffix}\".")

        self.status_msg = f"{status_msg_prefix}{self.suggested_model_part_name}{status_msg_suffix}"
        KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.GetRootModelPart(), self.status_msg)
        return f"{self.GetRootModelPart().FullName()}.{self.suggested_model_part_name}"

    def GetModelPart(self) -> Kratos.ModelPart:
        model_part_name = self.GetModelPartFullName()
        if self.status_msg == "" or self.model.HasModelPart(model_part_name):
            # if it is already there, that means it is already created. Hence return it or
            # if the self.status_msg == "" means, there is no need to do any model part operations.
            # therefore the model part should exist already.
            return self.model[model_part_name]
        else:
            # the model part with the required operation not found. Hence creating it.
            sub_model_part = self.model.CreateModelPart(model_part_name)
            # now fill the submodel part
            operation_model_parts = [self.model[name] for name in self.list_of_operation_model_part_full_names]
            if self.operation_type == ModelPartOperation.OperationType.UNION:
                Kratos.ModelPartOperationUtilities.Union(sub_model_part, self.root_model_part, operation_model_parts, self.add_neighbours)
            elif self.operation_type == ModelPartOperation.OperationType.INTERSECT:
                Kratos.ModelPartOperationUtilities.Intersect(sub_model_part, self.root_model_part, operation_model_parts, self.add_neighbours)

            Kratos.Logger.PrintInfo("ModelPartUtilities", f"Created sub model part \"{sub_model_part.FullName()}\".")

            return sub_model_part

class ModelPartUtilities:
    @staticmethod
    def __GenerateUniqueIdentifier(prefix: str, list_of_names: 'list[str]', add_neighbours: bool) -> str:
        if not all([name.find("#") == -1 for name in list_of_names]):
            raise RuntimeError(f"The provided model part names has \"#\" which is invalid. Please remove all occurances of \"#\" character. Provided model part names:\n\t" + "\n\t".join(list_of_names))

        sorted_names = sorted(list_of_names)
        post_fix = "IN" if add_neighbours else "EN"
        return (f"{prefix}_{post_fix}_" + "_".join(sorted_names)).replace(".", "-")

    @staticmethod
    def GetMergedMap(input_dict: 'dict[Any, KratosOA.CollectiveExpression]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            merging_model_part_names = [container_expression.GetModelPart().FullName() for container_expression in v.GetContainerExpressions()]

            if not merging_model_part_names:
                raise RuntimeError("Merging requires atleast one model part.")

            uniqe_identifier_name = ModelPartUtilities.__GenerateUniqueIdentifier("UNION", merging_model_part_names, add_neghbours)
            merged_model_part = ModelPartOperation(v.GetContainerExpressions()[0].GetModelPart().GetModel(), ModelPartOperation.OperationType.UNION, uniqe_identifier_name, merging_model_part_names, add_neghbours).GetModelPart()
            result[k] = merged_model_part

        return result

    @staticmethod
    def GetIntersectedMap(main_model_part: Kratos.ModelPart, input_dict: 'dict[Any, Kratos.ModelPart]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            intersecting_model_part_names = [main_model_part.FullName(), v.FullName()]
            uniqe_identifier_name = ModelPartUtilities.__GenerateUniqueIdentifier("INTERSECT", intersecting_model_part_names, add_neghbours)
            intersected_model_part = ModelPartOperation(main_model_part.GetModel(), ModelPartOperation.OperationType.INTERSECT, uniqe_identifier_name, intersecting_model_part_names, add_neghbours).GetModelPart()
            result[k] = intersected_model_part

        return result