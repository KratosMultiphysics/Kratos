import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
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
        self.FindIntersectionPointsAndSurfaces_Hexa3D8N(mp)

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


    def FindIntersectionPointsAndSurfaces_Hexa3D8N(self, mp_interface: Kratos.ModelPart) -> None:
        control_field = self.control_field.Evaluate()

        # Local edges of Hexahedra3D8N
        edge_pairs = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ]

        # Precompute node neighbours to detect boundary nodes
        find_neighbours = Kratos.FindGlobalNodalNeighboursProcess(self.model_part)
        find_neighbours.Execute()
        node_neighbours = find_neighbours.GetNeighbourIds(self.model_part.Nodes)
        expected_valence = 17  # internal node in 2D Hexa3D8 structured mesh

        element: Kratos.Element
        for element in self.model_part.Elements:
            geom = element.GetGeometry()
            zero_level_set = []
            corners = []

            # usual zero-crossing along edges
            for (i1, i2) in edge_pairs:
                n1 = geom[i1]
                n2 = geom[i2]
                phi1 = control_field[n1.Id - 1]
                phi2 = control_field[n2.Id - 1]

                # Sign change -> zero-crossing
                if phi1 * phi2 < 0.0:
                    p_int = self.InterpolateZeroCrossing(n1, n2, phi1, phi2)
                    zero_level_set.append(p_int)

            # include boundary nodes that are part of material
            for node in geom:
                phi = control_field[node.Id - 1]
                num_neigh = len(node_neighbours[node.Id])
                # print(f"NEIGHBOURS: {num_neigh}")
                if num_neigh == 11 and phi >= 0.0:
                    zero_level_set.append((node.X0, node.Y0, node.Z0))
                elif num_neigh == 7 and phi >= 0.0:
                    corners.append((node.X0, node.Y0, node.Z0))

            # create surfaces if any intersection points exist
            if len(zero_level_set) > 0:
                self.FindSurfaces(mp_interface, zero_level_set, corners)
    
    def InterpolateZeroCrossing(self, node1: Kratos.Node, node2: Kratos.Node, phi1: float, phi2: float):
        lmbd = phi1 / (phi1 - phi2)
        x = node1.X0 - lmbd * (node1.X0 - node2.X0)
        y = node1.Y0 - lmbd * (node1.Y0 - node2.Y0)
        z = node1.Z0 - lmbd * (node1.Z0 - node2.Z0)
        return (x, y, z)

    def FindSurfaces(self, mp_interface: Kratos.ModelPart, zero_level_set: list, corners: list = []) -> None:
        prop = mp_interface.GetProperties(2)
        i = mp_interface.NumberOfNodes()
        j = mp_interface.NumberOfConditions()
        n = len(zero_level_set)

        node_ids = []
        for k in range(n):
            i += 1
            mp_interface.CreateNewNode(i, zero_level_set[k][0], zero_level_set[k][1], zero_level_set[k][2])
            node_ids.append(i)
        if corners == []:
            if n == 2:
                mp_interface.CreateNewCondition("LineCondition3D2N", j + 1, node_ids, prop)
            elif n == 3:
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 1, node_ids, prop)
            elif n == 4:
                # Find triangulation
                distances = {}
                for l in range(3):
                    for k in range(l + 1, 4):
                        node1 = zero_level_set[l]
                        node2 = zero_level_set[k]
                        distance = sum((node1[dim] - node2[dim])**2 for dim in range(3))
                        distances[distance] = (node_ids[l], node_ids[k])
                furthest = max(distances)
                a, b = distances[furthest]

                # Collect all node indices
                all_nodes = set(node_ids)

                # Remaining two nodes (not in the longest diagonal)
                remaining = list(all_nodes - {a, b})
                c, d = remaining[0], remaining[1]

                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, c, b, d], prop)
                # # Create two triangle conditions sharing the shorter diagonal (c,d)
                # mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 1, [a, c, d], prop)
                # mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 2, [b, c, d], prop)
            else:
                Kratos.Logger.PrintWarning("InterfaceOutputProcess", f"Element produced {n} intersections — not handled.")
        else:
            if len(corners) != 2:
                raise RuntimeError("Expected exactly 2 corner nodes for corner elements.")
            corner_node_ids = []
            for k in range(2):
                i += 1
                mp_interface.CreateNewNode(i, corners[k][0], corners[k][1], corners[k][2])
                corner_node_ids.append(i)
            # corner case: split 6 nodes into 2 quads -> 4 triangles
            # raise RuntimeError(f"Essa: {coords}\n {node_ids}")
            distances = {}
            for k in range(1,4):
                node1 = 0
                node2 = k
                d = sum((zero_level_set[node1][dim] - zero_level_set[node2][dim]) ** 2 for dim in range(3))
                #save distance from node 0 to all other non corner nodes by id
                distances[d] = node_ids[node2]
            # 2 highest distances belong to another face, c is with checked node
            a = distances.pop(max(distances))
            b = distances.pop(max(distances))
            c = distances[max(distances)]
            d = node_ids[0]

            # Order each quad along longest diagonal
            mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, b, corner_node_ids[0], corner_node_ids[1]], prop)
            mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, [c, d, corner_node_ids[0], corner_node_ids[1]], prop)