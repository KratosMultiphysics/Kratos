import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from pathlib import Path
import csv

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
        self.output_path_cloud = "point_cloud"
        self.save_output_files_in_folder = parameters["save_output_files_in_folder"].GetBool()
        # self.output_interval = parameters["output_interval"].GetInt()
        self.model_part = model[parameters["model_part_name"].GetString()]
        self.optimization_problem = optimization_problem
        self.output_interval = parameters["output_interval"].GetInt()
        # self.created_nodes = {}
        self.created_coordinates = {}


    def IsOutputStep(self):
        if self.optimization_problem.GetStep() % self.output_interval == 0:
            return True
        else:
            return False

    def PrintOutput(self):
        if self.optimization_problem.GetStep() % self.output_interval == 0:
            for control in self.optimization_problem.GetListOfControls():
                self.control_field = control.GetControlField()
            model = Kratos.Model()
            mp = model.CreateModelPart("Interface")
            mp.CreateNewProperties(2)
            self.created_coordinates = {}
            self.FindIntersectionPointsAndSurfaces_Hexa3D8N(mp)

            if self.save_output_files_in_folder:
                self.output_path = Path(self.output_path)
                if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
                    kratos_utils.DeleteDirectoryIfExisting(str(self.output_path_cloud))
                    self.model_part.ProcessInfo[Kratos.IS_RESTARTED] = True
                self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
                # now create the output path
                Kratos.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path))
                Kratos.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path_cloud))
            else:
                self.output_path = Path(".")
                self.output_path_cloud = Path(".")

            output_file_name = self.output_file_name_prefix
            output_file_name = output_file_name.replace("<model_part_full_name>", self.model_part.FullName())
            output_file_name = output_file_name.replace("<model_part_name>", self.model_part.Name)
            output_file_name = output_file_name.replace("<step>", str(self.optimization_problem.GetStep()))
            vtu_output = Kratos.VtuOutput(mp)
            vtu_output.PrintOutput(str(self.output_path /output_file_name))
            self.write_created_coordinates_to_csv(self.output_path_cloud + "/interface_nodes_" + str(self.optimization_problem.GetStep()) + ".csv")


    def write_created_coordinates_to_csv(self, filename="created_coordinates.csv"):
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            # Write header
            writer.writerow(["node_id", "x", "y", "z", "type"])

            # Write each entry
            for (x, y, z), (node_id, node_type) in self.created_coordinates.items():
                writer.writerow([node_id, x, y, z, node_type])

        print(f"✅ CSV file written: {filename}")

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

        element: Kratos.Element
        for element in self.model_part.Elements:
            geom = element.GetGeometry()
            zero_level_set = []
            corners = []
            edge = []
            ntype = []

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
                    ntype.append("intersection")

            # include boundary nodes that are part of material
            #####
            cords = []
            ids = []
            ######
            for node in geom:
                #####
                cords.append((node.X0, node.Y0, node.Z0))
                ids.append(node.Id)
                ######
                phi = control_field[node.Id - 1]
                num_neigh = len(node_neighbours[node.Id])
                if num_neigh == 11 and phi >= 0.0:
                    ntype.append("edge")
                    zero_level_set.append((node.X0, node.Y0, node.Z0))
                # elif num_neigh == 17 and phi >= 0.0:
                #     ntype.append("surface")
                #     zero_level_set.append((node.X0, node.Y0, node.Z0))
                elif num_neigh == 7 and phi >= 0.0:
                    ntype.append("corner")
                    corners.append((node.X0, node.Y0, node.Z0))

            # create surfaces if any intersection points exist
            if len(zero_level_set) > 0:
                self.FindSurfaces(mp_interface, zero_level_set, corners, edge, cords, ids, ntype)
    
    def InterpolateZeroCrossing(self, node1: Kratos.Node, node2: Kratos.Node, phi1: float, phi2: float):
        lmbd = phi1 / (phi1 - phi2)
        x = node1.X0 - lmbd * (node1.X0 - node2.X0)
        y = node1.Y0 - lmbd * (node1.Y0 - node2.Y0)
        z = node1.Z0 - lmbd * (node1.Z0 - node2.Z0)
        return (x, y, z)

    def FindSurfaces(self, mp_interface: Kratos.ModelPart, zero_level_set: list, corners: list = [], edge: list = [], temp1: list = [], temp2: list = [], ntype = []) -> None:
        prop = mp_interface.GetProperties(2)
        i = mp_interface.NumberOfNodes()
        j = mp_interface.NumberOfConditions()
        n = len(zero_level_set)

        node_ids = []
        for k in range(n):
            if not zero_level_set[k] in self.created_coordinates:
                i += 1
                mp_interface.CreateNewNode(i, zero_level_set[k][0], zero_level_set[k][1], zero_level_set[k][2])
                node_ids.append(i)
                self.created_coordinates[zero_level_set[k]] = (i, ntype[k])
                # self.created_nodes[i] = zero_level_set[k]
            else:
                id = self.created_coordinates[zero_level_set[k]][0]
                node_ids.append(id)

        # print(f"id:{node_ids}\ncoordinates: {zero_level_set}")

        if corners == []:
            # raise RuntimeError(f"zero node coordinates: {zero_level_set}\nzero node ids: {node_ids}\ncorners: {corners}\nedge: {edge}")
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
            # 5 node internal intersection
            elif n == 5:
                node1 = zero_level_set[0]
                node1_id = node_ids[0]
                surface = [node1_id]
                for k in range(len(zero_level_set)-1):
                    del zero_level_set[node_ids.index(node1_id)]
                    node_ids.remove(node1_id)
                    node1_id = self.FindClosest(zero_level_set, node_ids, node1)[0]
                    node1 = zero_level_set[node_ids.index(node1_id)]
                    surface.append(node1_id)
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 1, surface[0:3], prop)
                surface.pop(1)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, surface, prop)
            elif n == 6 and len(edge)!=0:
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
            elif n == 6 and len(edge) == 0:
                shortest = []
                third_shortest = {}
                for l in range(6):
                    distances = {}
                    node1 = zero_level_set[l]
                    for k in range(6):
                        if l == k:
                            continue
                        node2 = zero_level_set[k]
                        distance = sum((node1[dim] - node2[dim])**2 for dim in range(3))
                        distances[(node_ids[l], node_ids[k])] = distance
                    cur_shortest = sorted(distances, key=distances.get, reverse=False)[:3]
                    third_shortest[cur_shortest[2]] = distances[cur_shortest[2]]
                    if (cur_shortest[0] not in shortest) and (cur_shortest[0][::-1] not in shortest):
                        shortest.append(cur_shortest[0])
                    if cur_shortest[1] not in shortest and cur_shortest[1][::-1] not in shortest:
                        shortest.append(cur_shortest[1])
                if len(shortest) == 6:
                    third = sorted(third_shortest, key=third_shortest.get, reverse=False)[0]
                    shortest.append(third)
                # how many instances
                shared_nodes = {}
                for pair in shortest:
                    for node in pair:
                        if not node in shared_nodes:
                            shared_nodes[node] = 1
                        else:
                            shared_nodes[node] = shared_nodes[node] + 1
                shared = []
                for shared_node, instances in shared_nodes.items():
                    if instances == 3:
                        shared.append(shared_node)
                if tuple(shared) in shortest:
                    shortest.remove(tuple(shared))
                elif tuple(shared[::-1]) in shortest:
                    shortest.remove(tuple(shared[::-1]))
                surf1 = shared.copy()
                cur = shared[1]
                while len(surf1) < 4:
                    for l in range(len(shortest)):
                        if cur == shortest[l][0]:
                            surf1.append(shortest[l][1])
                            cur = shortest[l][1]
                            shortest.pop(l)
                            break
                        elif cur == shortest[l][1]:
                            surf1.append(shortest[l][0])
                            cur = shortest[l][0]
                            shortest.pop(l)
                            break
                if (surf1[0], surf1[-1]) in shortest:
                    shortest.remove((surf1[0], surf1[-1]))
                elif (surf1[-1], surf1[0]) in shortest:
                    shortest.remove((surf1[-1], surf1[0]))
                surf2 = shared.copy()
                cur = shared[1]
                while len(surf2) < 4:
                    for l in range(len(shortest)):
                        if cur == shortest[l][0]:
                            surf2.append(shortest[l][1])
                            cur = shortest[l][1]
                            shortest.pop(l)
                            break
                        elif cur == shortest[l][1]:
                            surf2.append(shortest[l][0])
                            cur = shortest[l][0]
                            shortest.pop(l)
                            break
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, surf1, prop)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, surf2, prop)
            elif n == 7:
                Kratos.Logger.PrintWarning("InterfaceOutputProcess", f"Element produced {n} intersections — not handled. Hiljem lisan")
            # two internal node intersections
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
            # intersection and boundary element
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
                print(f"element Coordinates: {temp1}\n Element ids: {temp2}")
                
                raise RuntimeError(f"Node ids: {node_ids}\nCoordinates: {zero_level_set}\n Edges: {edge}")
        else:
            n = len(corners)
            corner_node_ids = []
            for k in range(n):
                i += 1
                mp_interface.CreateNewNode(i, corners[k][0], corners[k][1], corners[k][2])
                corner_node_ids.append(i)
            # triangular cut including corner (??)
            if n + len(zero_level_set) == 3:
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", j + 1, node_ids + corner_node_ids, prop)
            # corner case: split 6 nodes into 2 quads -> 4 triangles
            # raise RuntimeError(f"Essa: {coords}\n {node_ids}")
            # pseudo 3D full corner
            elif len(zero_level_set) == 4 and n == 2:
                res = self.FindClosest(zero_level_set, node_ids, zero_level_set[0])
                # print(f"corner nodes: {corner_node_ids}, Zero-level: {distances}")
                a, b, c = res[2], res[1], res[0]
                d = node_ids[0]
                # Order each quad along longest diagonal
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 1, [a, b, corner_node_ids[0], corner_node_ids[1]], prop)
                mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + 2, [c, d, corner_node_ids[0], corner_node_ids[1]], prop)
            # 3D full corner
            elif len(zero_level_set) == 6 and n == 1:
                corner_closest = self.FindClosest(zero_level_set, node_ids, corners[0])[:3]
                for k, lvl_one in enumerate(corner_closest):
                    lvl_two = self.FindClosest(zero_level_set, node_ids, zero_level_set[node_ids.index(lvl_one)])[0]
                    lvl_three = self.FindClosest(zero_level_set, node_ids, zero_level_set[node_ids.index(lvl_two)])
                    if lvl_three[0] == lvl_one:
                        lvl_three = lvl_three[1]
                    else:
                        lvl_three = lvl_three[0]
                    mp_interface.CreateNewCondition("SurfaceCondition3D4N", j + k + 1, [corner_node_ids[0], lvl_one, lvl_two, lvl_three], prop)
                    del zero_level_set[node_ids.index(lvl_two)]
                    node_ids.remove(lvl_two)
            else:
                print(f"corner nodes: {corner_node_ids}, Zero-level: {node_ids}")

    def FindClosest(self, searced_node_list, node_ids, anchor_node):
        distances = {}
        for k in range(0, len(searced_node_list)):
            node2 = k
            if anchor_node == searced_node_list[k]:
                continue
            d = sum((anchor_node[dim] - searced_node_list[node2][dim]) ** 2 for dim in range(3))
            #save distance from node 0 to all other non corner nodes by id
            distances[node_ids[node2]] = d
        # return sorted ids of closest nodes
        return sorted(distances, key=distances.get, reverse=False)