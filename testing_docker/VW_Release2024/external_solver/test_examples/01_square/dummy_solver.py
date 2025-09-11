from kratos_external_solver_optimization import external_solver
from kratos_external_solver_optimization import model
from kratos_external_solver_optimization import mesh_writer
import json
from distutils.dir_util import copy_tree
import os
import numpy as np

def write_mdpa(mdpa_file_to_be_created, kratos_model):
    writer = mesh_writer.MeshWriter(mdpa_file_to_be_created, kratos_model)

    with open(mdpa_file_to_be_created, 'w') as out_file:
        writer._write_properties(out_file)
        writer._write_nodes(out_file)
        writer._write_elements(out_file)
        writer._write_conditions(out_file)
        writer._write_sub_model_parts(out_file)

class DummySolver(external_solver.ExternalSolver):

    def __init__(self, parameter_file_name):
        self.solver_name = "dummy_solver"
        parameter_file = open(parameter_file_name)
        self.parameter = json.load(parameter_file)
        parameter_file.close()

        self.dummy_model = None

        # parameter für quadratische geometrie
        self.number_nodes_x = self.parameter["solver_settings"]["number_of_nodes_x"]
        self.number_nodes_y = self.parameter["solver_settings"]["number_of_nodes_y"]

    def GetCaseName(self) -> str:
        return self.parameter["solver_settings"]["case_name"]

    def GetOriginalFilePath(self) -> str:
        return self.parameter["solver_settings"]["file_path"]

    def InitializeBeforeOptimization(self) -> None:
        """
        Erstelle quadratische Geometrie als Dummy Model
        und speichere diese als .mdpa Datei im Ordner
        der originalen Solver Dateien ('/Dummy/')
        """
        NodeDict, ElementDict = self._create_quadratic_mesh(self.number_nodes_x, self.number_nodes_y)
        self.dummy_model = model.Model("test")
        for id, coordinates in NodeDict.items():
            self.dummy_model.create_node(id, coordinates[0], coordinates[1], coordinates[2])

        self.dummy_model.create_material("1", "1", "1", "1")
        self.dummy_model.create_properties("1", "1", "dummy_property", "1")

        for id, node_ids in ElementDict.items():
            self.dummy_model.create_condition(id, 1, node_ids, "SurfaceCondition3D4N")

        self.dummy_model.add_sub_model_part("design_surface")
        design_sub_model = self.dummy_model.sub_model_parts["design_surface"]

        for id in NodeDict.keys():
            design_sub_model.add_node(id)

        for id in ElementDict.keys():
            design_sub_model.add_condition(id)

        edge_ids = []
        for id, coordinates in NodeDict.items():
            if coordinates[0] == 1 or coordinates[0] == self.number_nodes_x:
                edge_ids.append(id)
                continue
            if coordinates[1] == 0 or coordinates[1] == self.number_nodes_y - 1:
                edge_ids.append(id)
                continue

        self.dummy_model.add_sub_model_part("edges")
        edge_sub_model = self.dummy_model.sub_model_parts["edges"]
        for id in edge_ids:
            edge_sub_model.add_node(id)

        write_mdpa(f"Dummy/{self.GetCaseName()}.mdpa", self.dummy_model)

    def TranslateToKratosModel(self) -> None:
        """
        Übersetze Dummy Model in .mdpa Datei für die
        Kratos Formoptimierung und speichere diese
        im Hauptordner ab.
        """
        write_mdpa("structure.mdpa", self.dummy_model)

    def CopyInputFilesToCurrentIteration(self, current_iteration_path,
                                         previous_iteration_path,
                                         original_path,
                                         iteration) -> None:
        """
        Kopiere Solver Dateien in derzeitigen Optimierungspfad.
        Im ersten Optimierungsschritt werden die Dateien aus dem
        originalen Ordern ('/Dummy/') kopiert.
        """
        if iteration == 1:
            copy_tree(original_path, current_iteration_path)
            return

        copy_tree(previous_iteration_path, current_iteration_path)

    def UpdateMesh(self, current_iteration_path,
                   previous_iteration_path,
                   iteration,
                   original_path,
                   shape_update, thickness_update,
                   shape, thickness, current_design) -> None:
        """
        Aktualisiere Dummy Model und überschreibe die
        .mdpa Datei des Dummy Solvers im derzeitigen
        Optimierungspfad mit den neuen Informationen.
        """
        design_sub_model = self.dummy_model.sub_model_parts["design_surface"]
        for node in design_sub_model.nodes.values():
            node.x += shape_update[node.id][0]
            node.y += shape_update[node.id][1]
            node.z += shape_update[node.id][2]

        mdpa_file_path =  os.path.join(current_iteration_path, f"{self.GetCaseName()}.mdpa")

        write_mdpa(mdpa_file_path, self.dummy_model)

    def RunAnalysis(self,
                    current_iteration_path,
                    current_design,
                    iteration) -> None:
        pass

    def ReadValue(self,
                  identifier,
                  current_iteration_path,
                  current_design,
                  iteration) -> float:
        return 100.0 - iteration

    def ReadShapeGradient(self,
                          identifier,
                          current_iteration_path,
                          current_design,
                          iteration) -> dict:
        """
        Nur der Mittelknoten erhält eine Sensitivität
        in negativer z-Richtung.
        """
        center_node_id = int(np.ceil(self.number_nodes_x * self.number_nodes_y / 2))
        shape_gradient = dict()
        for node in self.dummy_model.nodes.values():
            gradient = [0.0, 0.0, 0.0]
            if node.id == center_node_id:
                gradient = [0.0, 0.0, -1.0]
            shape_gradient[node.id] = gradient

        return shape_gradient

    def ReadThicknessGradient(self,
                                identifier,
                                current_iteration_path,
                                current_design,
                                iteration) -> dict:
        pass

    def _create_quadratic_mesh(self, nNodesX, nNodesY):
        """
        Erstellt uniformes quadratisches Netz.
        """
        NodeDict = {}
        yRow = 0
        for i in range(1, nNodesX * nNodesY + 1, 1):
            xColumn = i - nNodesX * yRow
            print(f"xColumn: {xColumn}")
            if xColumn > nNodesX:
                yRow += 1
                xColumn = i - nNodesX * yRow
            tmpx = xColumn
            tmpy = yRow
            tmpz = 0.0

            print(f"tmp: [{tmpx}, {tmpy}, {tmpz}]")
            NodeDict[i] = [tmpx, tmpy, tmpz]

        ElementDict = {}
        yRow = 0
        nElementsX = nNodesX - 1
        nElementsY = nNodesY - 1
        for i in range(1, nElementsX * nElementsY + 1):
            xColumn = i - nElementsX * yRow
            if xColumn > nElementsX:
                yRow += 1
                xColumn = i - nElementsX * yRow

            lo = (yRow) * nNodesX + xColumn
            ro = lo + 1
            lu = lo + nNodesX
            ru = lu + 1

            ElementDict[i] = [lo, lu, ru, ro]

        return NodeDict, ElementDict
