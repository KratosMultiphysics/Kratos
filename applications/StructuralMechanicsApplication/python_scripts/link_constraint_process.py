# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication

# --- STD Imports ---
import typing


class LinkConstraintProcess(KratosMultiphysics.Process):
    """@brief Create and manage @ref LinkConstraint "constraints" between pairs of @ref Node "nodes" that force them to keep the distance between each other constant."""


    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__model_part: KratosMultiphysics.ModelPart = model.GetModelPart(parameters["model_part_name"].GetString())
        self.__node_pairs: "list[tuple[int,int]]" = []
        self.__dimensions: int = parameters["dimensions"].GetInt()
        self.__is_mesh_moved: bool = parameters["move_mesh_flag"].GetBool()
        self.__interval: KratosMultiphysics.IntervalUtility = KratosMultiphysics.IntervalUtility(parameters)
        self.__is_active: bool = False
        self.__constraints: "list[KratosMultiphysics.StructuralMechanicsApplication.LinkConstraint]" = []

        pair: KratosMultiphysics.Parameters
        for pair in parameters["node_pairs"].values():
            self.__node_pairs.append((pair[0].GetInt(), pair[1].GetInt()))


    @classmethod
    def GetDefaultParameters(cls: "typing.Type[LinkConstraintProcess]") -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters(R"""{
            "model_part_name" : "",
            "node_pairs" : [],
            "dimensions" : 3,
            "move_mesh_flag" : false,
            "interval" : [0, "End"]
        }""")


    def ExecuteBeforeSolutionLoop(self) -> None:
        self.__is_active = self.__interval.IsInInterval(self.__model_part.ProcessInfo[KratosMultiphysics.TIME])

        id: int = self.__GetLastConstraintId() + 1
        for id_first, id_second in self.__node_pairs:
            constraint = KratosMultiphysics.StructuralMechanicsApplication.LinkConstraint(
                id,
                self.__model_part.GetNode(id_first),
                self.__model_part.GetNode(id_second),
                self.__dimensions,
                self.__is_mesh_moved)

            constraint.Set(KratosMultiphysics.ACTIVE, self.__is_active)

            self.__constraints.append(constraint)
            self.__model_part.AddMasterSlaveConstraint(constraint)
            id += 1


    def ExecuteInitializeSolutionStep(self) -> None:
        is_active: bool = self.__interval.IsInInterval(self.__model_part.ProcessInfo[KratosMultiphysics.TIME])
        if is_active != self.__is_active:
            self.__is_active = is_active
            for constraint in self.__constraints:
                constraint.Set(KratosMultiphysics.ACTIVE, self.__is_active)


    def __GetLastConstraintId(self) -> int:
        constraint_count: int = len(self.__model_part.MasterSlaveConstraints)
        last_constraint_id: int
        if constraint_count:
            last_constraint_id = self.__model_part.MasterSlaveConstraints.back().Id
        else:
            last_constraint_id = 0
        return self.__model_part.GetCommunicator().GetDataCommunicator().MaxAll(last_constraint_id)


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> LinkConstraintProcess:
    return LinkConstraintProcess(model, parameters["Parameters"])
