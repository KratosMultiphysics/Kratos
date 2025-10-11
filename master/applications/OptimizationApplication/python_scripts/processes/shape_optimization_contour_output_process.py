import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics import FindGlobalNodalNeighboursProcess
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from pathlib import Path

def Factory(Model: Kratos.Model, parameters: Kratos.Parameters,  optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if(type(parameters) != Kratos.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return InterfaceOutputProcess(Model, parameters["settings"], optimization_problem)

class InterfaceOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>_<step>",
                "output_path"                 : "Solid_interface",
                "save_output_files_in_folder" : true,
                "output_interval"             : 1,
                "model_part_name"             : ""
            }
            """
        )
    
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.file_name = parameters["file_name"].GetString()
        self.output_file_name_prefix = parameters["file_name"].GetString()
        self.output_path = parameters["output_path"].GetString()
        self.save_output_files_in_folder = parameters["save_output_files_in_folder"].GetBool()
        # self.output_interval = parameters["output_interval"].GetInt()
        self.model_part = model[parameters["model_part_name"].GetString()]
        self.optimization_problem = optimization_problem


    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        for control in self.optimization_problem.GetListOfControls():
            self.control_field = control.GetControlField()
        model = Kratos.Model()
        mp = model.CreateModelPart("Interface")
        mp.CreateNewProperties(2)
        self.FindIntersectionPointsAndSurfaces2(mp)
        # mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        # mp.CreateNewNode(2, 1.0, 1.0, 0.0)
        # mp.CreateNewNode(3, 0.5, 0.5, 0.0)

        # mp.CreateNewCondition("LineCondition3D2N", 1, [1, 2], prop)
        # mp.CreateNewCondition("LineCondition3D2N", 2, [2, 3], prop)

        if self.save_output_files_in_folder:
            self.output_path = Path(self.output_path)
            if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
                self.model_part.ProcessInfo[Kratos.IS_RESTARTED] = True
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
            # now create the output path
            Kratos.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path))
        else:
            self.output_path = Path(".")

        output_file_name = self.output_file_name_prefix
        output_file_name = output_file_name.replace("<model_part_full_name>", self.model_part.FullName())
        output_file_name = output_file_name.replace("<model_part_name>", self.model_part.Name)
        output_file_name = output_file_name.replace("<step>", str(self.optimization_problem.GetStep()))
        vtu_output = Kratos.VtuOutput(mp)
        vtu_output.PrintOutput(str(self.output_path /output_file_name))

# with neighbours but gives more nodes than I expected
    def FindIntersectionPointsAndSurfaces(self, mp_interface: Kratos.ModelPart) -> None:
        control_field = Kratos.Expression.NodalExpression(self.model_part)
        control_field = self.control_field.Evaluate()

        # get a dict of neighbouring nodes
        find_neighbours = FindGlobalNodalNeighboursProcess(self.model_part)
        find_neighbours.Execute()
        node_neighbours = find_neighbours.GetNeighbourIds(self.model_part.Nodes)

        element: Kratos.Element
        for element in self.model_part.Elements:
            zero_level_set = []
            node: Kratos.Node
            element_nodes = element.GetNodes()
            for node in element.GetNodes():
                # check if it is already a point of level set. 
                # Add an error delta please.
                if control_field[node.Id - 1] == 0:
                    intersection_point = (node.X0, node.Y0, node.Z0)
                    zero_level_set.append(intersection_point)
                    # remove from list so we don't check again
                    element_nodes.remove(node)
                    continue
                # check line segments of current node
                neighbours = node_neighbours[node.Id]
                for node2 in element_nodes:
                    print(f"Current element node: {node.Id}, Checked node: {node2.Id}, Neighbours of current: {neighbours}")
                    if node2.Id in neighbours:
                        print("In neighbours")
                        if control_field[node.Id - 1] * control_field[node2.Id - 1] < 0:
                            zero_level_set.append(self.FindIntersectionPoints(node, node2))
                    else:
                        print("Not in neighbours")
                # remove from list so we don't check again 
                # (we already checked all of the line segments in this element that include current node)
                element_nodes.remove(node)
                # go to next node of this element
            # now we have a list of level set points of current element and can create a surface inbetween the points
            # print(f"Intersection points: {len(zero_level_set)}")
            if len(zero_level_set) > 0:
                self.FindSurfaces(mp_interface, zero_level_set)

    def FindIntersectionPointsAndSurfaces2(self, mp_interface: Kratos.ModelPart) -> None:
        control_field = Kratos.Expression.NodalExpression(self.model_part)
        control_field = self.control_field.Evaluate()

        # get a dict of neighbouring nodes
        find_neighbours = FindGlobalNodalNeighboursProcess(self.model_part)
        find_neighbours.Execute()
        node_neighbours = find_neighbours.GetNeighbourIds(self.model_part.Nodes)

        element: Kratos.Element
        for element in self.model_part.Elements:
            zero_level_set = []
            node: Kratos.Node
            element_nodes = element.GetNodes()
            for node in element.GetNodes():
                # check if it is already a point of level set. 
                # Add an error delta please.
                if control_field[node.Id - 1] == 0:
                    intersection_point = (node.X0, node.Y0, node.Z0)
                    zero_level_set.append(intersection_point)
                    # remove from list so we don't check again
                    element_nodes.remove(node)
                    continue
                # check line segments of current node
                for node2 in element_nodes:
                    if ((node.X0 == node2.X0 and node.Y0 == node2.Y0) or \
                        (node.X0 == node2.X0 and node.Z0 == node2.Z0) or \
                        (node.Z0 == node2.Z0 and node.Y0 == node2.Y0)) and \
                        (node != node2):
                        if control_field[node.Id - 1] * control_field[node2.Id - 1] < 0:
                            zero_level_set.append(self.FindIntersectionPoints(node, node2))
                # remove from list so we don't check again 
                # (we already checked all of the line segments in this element that include current node)
                element_nodes.remove(node)
                # go to next node of this element
            # now we have a list of level set points of current element and can create a surface inbetween the points
            # print(f"Intersection points: {len(zero_level_set)}")
            if len(zero_level_set) > 0:
                self.FindSurfaces(mp_interface, zero_level_set)

    def FindIntersectionPoints(self, node1: Kratos.Node, node2: Kratos.Node) -> tuple:
        control_field = Kratos.Expression.NodalExpression(self.model_part)
        control_field = self.control_field.Evaluate()
        d_phi = abs(control_field[node1.Id - 1] - control_field[node2.Id - 1])
        d_x = abs(node1.X0 - node2.X0) / d_phi * abs(control_field[node1.Id - 1])
        d_y = abs(node1.Y0 - node2.Y0) / d_phi * abs(control_field[node1.Id - 1])
        d_z = abs(node1.Z0 - node2.Z0) / d_phi * abs(control_field[node1.Id - 1])
        intersection_point = (node1.X0 + d_x, node1.Y0 + d_y, node1.Z0 + d_z)
        return intersection_point

    def FindSurfaces(self, mp_interface: Kratos.ModelPart, zero_level_set: list) -> None:
        prop = mp_interface.GetProperties(2)
        i = mp_interface.NumberOfNodes()
        j = mp_interface.NumberOfConditions()

        if len(zero_level_set) == 8:
            # all values are 0, first step 
            return 0

        for k in range(len(zero_level_set)):
            mp_interface.CreateNewNode(i + k + 1, zero_level_set[k][0], zero_level_set[k][1], zero_level_set[k][2])
            print(mp_interface.GetNode(i+k+1))

        if len(zero_level_set) == 2:
            mp_interface.CreateNewCondition("LineCondition3D2N", j + 1, [i + 1, i + 2], prop)
        elif len(zero_level_set) == 3:
            mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 1, [i + 1, i + 2, i + 3], prop)
        elif len(zero_level_set) == 4:
            mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [i + 1, i + 2, i + 3, i + 4], prop)
        else:
            raise RuntimeError(f"More than 4 intersection points in one element ({len(zero_level_set)} points). No implementation yet.")