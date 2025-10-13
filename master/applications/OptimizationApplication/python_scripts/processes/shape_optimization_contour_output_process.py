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
            edge = []

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
                    edge.append([(i1,i2), (phi1, phi2)])

            # include boundary nodes that are part of material
            for node in geom:
                phi = control_field[node.Id - 1]
                num_neigh = len(node_neighbours[node.Id])
                if num_neigh == 11 and phi >= 0.0:
                    zero_level_set.append((node.X0, node.Y0, node.Z0))
                elif num_neigh == 7 and phi >= 0.0:
                    corners.append((node.X0, node.Y0, node.Z0))

            # create surfaces if any intersection points exist
            if len(zero_level_set) > 0:
                self.FindSurfaces(mp_interface, zero_level_set, corners, edge)
    
    def InterpolateZeroCrossing(self, node1: Kratos.Node, node2: Kratos.Node, phi1: float, phi2: float):
        lmbd = phi1 / (phi1 - phi2)
        x = node1.X0 - lmbd * (node1.X0 - node2.X0)
        y = node1.Y0 - lmbd * (node1.Y0 - node2.Y0)
        z = node1.Z0 - lmbd * (node1.Z0 - node2.Z0)
        return (x, y, z)

    def FindSurfaces(self, mp_interface: Kratos.ModelPart, zero_level_set: list, corners: list = [], edge: list = []) -> None:
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
                        distances[(node_ids[l], node_ids[k])] = distance
                furthest = max(distances, key=distances.get)
                a, b = furthest

                # Remaining two nodes (not in the longest diagonal)
                remaining = list(set(node_ids) - {a, b})
                c, d = remaining[0], remaining[1]

                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, c, b, d], prop)
                # # Create two triangle conditions sharing the shorter diagonal (c,d)
                # mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 1, [a, c, d], prop)
                # mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 2, [b, c, d], prop)
            elif n == 6:
                distances = {}
                distances_bound = {}
                for l in range(4):
                    node1 = zero_level_set[l]
                    for k, node_bound in enumerate(zero_level_set[4:6]):
                        dist_bound = sum((node1[dim] - node_bound[dim])**2 for dim in range(3))
                        distances_bound[(node_ids[l], node_ids[k + 4])] = dist_bound
                    for k in range(l + 1, 4):
                        node2 = zero_level_set[k]
                        distance = sum((node1[dim] - node2[dim])**2 for dim in range(3))
                        distances[(node_ids[l], node_ids[k])] = distance
                        
                furthest = max(distances, key=distances.get)
                a, b = furthest
                # Collect all cutting node indices
                remaining = list(set(node_ids[0:4]) - {a, b})
                c, d = remaining[0], remaining[1]
                # boundary nodes
                closest = min(distances_bound, key=distances_bound.get)
                e, f = closest
                distances_bound.pop(closest)
                g, h = min(distances_bound, key=distances_bound.get)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, c, b, d], prop)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, [e, f, h, g], prop)
            elif n == 8 and len(edge) == 8:
                connections = []
                phi = 0
                for l in range(7):
                    for k in range(l + 1, 8):
                        common = set(edge[l][0]).intersection(edge[k][0])
                        if common:
                            if list(common)[0] == edge[l][0][0]:
                                phi = edge[l][1][0]
                            elif list(common)[0] == edge[l][0][1]:
                                phi = edge[l][1][1]

                            if phi > 0:
                                connections.append((l, k))
                        # Find pairs in z
                        if abs(zero_level_set[l][0] - zero_level_set[k][0]) < 1e-10 and \
                           abs(zero_level_set[l][1] - zero_level_set[k][1]) < 1e-10:
                            connections.append((l, k))
                
                surf1 = [connections[0][0], connections[0][1]]
                cur = connections[0][1]
                connections.pop(0)
                while len(surf1) < 4:
                    for l in range(len(connections)):
                        if cur == connections[l][0]:
                            surf1.append(connections[l][1])
                            cur = connections[l][1]
                            connections.pop(l)
                            break
                        elif cur == connections[l][1]:
                            surf1.append(connections[l][0])
                            cur = connections[l][0]
                            connections.pop(l)
                            break
                connections.remove((surf1[0], surf1[-1]))

                surf2 = [connections[0][0], connections[0][1]]
                cur = connections[0][1]
                connections.pop(0)
                while len(surf2) < 4:
                    for l in range(len(connections)):
                        if cur == connections[l][0]:
                            surf2.append(connections[l][1])
                            cur = connections[l][1]
                            connections.pop(l)
                            break
                        elif cur == connections[l][1]:
                            surf2.append(connections[l][0])
                            cur = connections[l][0]
                            connections.pop(l)
                            break
                            
                a, b, c, d = node_ids[surf1[0]], node_ids[surf1[1]], node_ids[surf1[2]], node_ids[surf1[3]]
                e, f, g, h = node_ids[surf2[0]], node_ids[surf2[1]], node_ids[surf2[2]], node_ids[surf2[3]]
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, b, c, d], prop)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, [e, f, g, h], prop)
            elif n == 8 and len(edge) == 4:
                distances = {}
                distances2 = {}
                for l in range(3):
                    for k in range(l + 1, 4):
                        node1 = zero_level_set[l]
                        node2 = zero_level_set[k]
                        node3 = zero_level_set[l+4]
                        node4 = zero_level_set[k+4]
                        distance = sum((node1[dim] - node2[dim])**2 for dim in range(3))
                        distances[(node_ids[l], node_ids[k])] = distance
                        distance2 = sum((node3[dim] - node4[dim])**2 for dim in range(3))
                        distances2[(node_ids[l+4], node_ids[k+4])] = distance2
                furthest = max(distances, key=distances.get)
                a, b = furthest
                furthest2 = max(distances2, key=distances2.get)
                e, f = furthest2
                # Remaining two nodes (not in the longest diagonal)
                remaining = list(set(node_ids[0:4]) - {a, b})
                c, d = remaining[0], remaining[1]
                remaining = list(set(node_ids[4:8]) - {e, f})
                g, h = remaining[0], remaining[1]

                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, c, b, d], prop)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, [e, g, f, h], prop)
            else:
                Kratos.Logger.PrintWarning("InterfaceOutputProcess", f"Element produced {n} intersections — not handled.")
                control_field = self.control_field.Evaluate()
                phi = []
                for id in node_ids:
                    phi.append(control_field[id])
                raise RuntimeError(f"Node ids: {node_ids}\nCoordinates: {zero_level_set}\nPhi value: {phi}")
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
                distances[node_ids[node2]] = d
            # 2 highest distances belong to another face, c is with checked node
            res = sorted(distances, key=distances.get, reverse=True)
            a, b, c = res[0], res[1], res[2]
            d = node_ids[0]

            # Order each quad along longest diagonal
            mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, b, corner_node_ids[0], corner_node_ids[1]], prop)
            mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, [c, d, corner_node_ids[0], corner_node_ids[1]], prop)