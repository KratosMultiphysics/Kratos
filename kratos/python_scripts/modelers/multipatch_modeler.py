#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Construct independent P1/Q1 meshes and interfaces, without physics."""

import math

import KratosMultiphysics as KM


def Factory(model, settings):
    return FemMultipatchModeler(model, settings)


class FemMultipatchModeler(KM.Modeler):
    """Build independent fitted boxes or conservative Q1 Gap-SBM interfaces.

    Meshes are submodelparts with globally distinct node and element IDs.
    Interfaces live outside the mesh submodelparts; gap mode also stores
    integration geometries carrying existing background nodes. True contours have their own nodes; surrogate edges
    reference background nodes. No conditions, DOFs or constraints are added.
    """

    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(KM.Parameters('''{
            "model_part_name": "FemMultipatch",
            "interface_mode": "fitted",
            "background": {},
            "internal_patches": [],
            "tolerance": 1e-10
        }'''))
        self.tolerance = self.settings["tolerance"].GetDouble()
        if not math.isfinite(self.tolerance) or self.tolerance <= 0.0:
            raise ValueError("tolerance must be finite and positive.")

        self.interface_mode = self.settings["interface_mode"].GetString()
        if self.interface_mode not in ("fitted", "gap_sbm"):
            raise ValueError("interface_mode must be 'fitted' or 'gap_sbm'.")

        # Validate box settings before creating entities. Gap ownership is
        # checked in C++ against the generated meshes before committing flags.
        background = self._read_box(self.settings["background"])
        if abs(background["rotation_angle"].GetDouble()) > self.tolerance:
            raise ValueError("The background must be Cartesian and unrotated.")
        if self.interface_mode == "gap_sbm" and background["element_type"].GetString() != "quadrilateral":
            raise ValueError("Gap-SBM initially requires Q1 quadrilaterals.")
        self.patch_settings = [background]
        self.boxes = [self._bounds(background)]
        for patch in self.settings["internal_patches"].values():
            patch = patch.Clone()
            if patch.Has('cylinder'): # to delete
                cylinder = patch['cylinder']
                center = list(cylinder['center'].GetVector())
                radius = cylinder['radius'].GetDouble()
                clearance = cylinder['clearance'].GetDouble()
                dimensions = list(patch['dimensions'].GetVector())
                if (len(center) != 2 or len(dimensions) != 2 or
                        not all(math.isfinite(x) for x in center + dimensions + [radius, clearance]) or
                        radius <= 0.0 or clearance < 0.0 or
                        min(dimensions) <= 2.0 * (radius + clearance)):
                    raise ValueError('Cylinder and clearance must fit strictly inside the inner patch.')
                origin = KM.Vector([center[i] - 0.5 * dimensions[i] for i in range(2)])
                if patch.Has('origin'):
                    patch['origin'].SetVector(origin)
                else:
                    patch.AddVector('origin', origin)
                patch.RemoveValue('cylinder')
            patch = self._read_box(patch)
            if self.interface_mode == "fitted":
                if abs(patch["rotation_angle"].GetDouble()) > self.tolerance:
                    raise ValueError("Rotation requires interface_mode='gap_sbm'.")
                self._align_internal_box(patch, background)
                bounds = self._bounds(patch)
                self._check_internal_box(bounds)
            else:
                if patch["element_type"].GetString() != "quadrilateral":
                    raise ValueError("Gap-SBM initially requires Q1 quadrilaterals.")
                width, height = patch["dimensions"].GetVector()
                if abs(width - height) > self.tolerance:
                    raise ValueError("The initial Gap-SBM inner patch must be square.")
                corners = self._rotated_corners(patch)
                bounds = (min(p[0] for p in corners), min(p[1] for p in corners),
                          max(p[0] for p in corners), max(p[1] for p in corners))
                bx0, by0, bx1, by1 = self.boxes[0]
                if not (bx0 + self.tolerance < bounds[0] and bounds[2] < bx1 - self.tolerance and
                        by0 + self.tolerance < bounds[1] and bounds[3] < by1 - self.tolerance):
                    raise ValueError("Rotated patches must lie strictly inside the background.")
            self.patch_settings.append(patch)
            self.boxes.append(bounds)

        name = self.settings["model_part_name"].GetString()
        if not name or "." in name:
            raise ValueError(
                "model_part_name must be a nonempty root model part name.")
        if model.HasModelPart(name):
            raise ValueError(f"Model part '{name}' already exists.")
        self.model_part = model.CreateModelPart(name)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 2
        self.patches = []

    def SetupModelPart(self):
        """Generate meshes, mark holes and store both boundaries."""
        super().SetupModelPart()
        if self.patches:
            raise RuntimeError(
                "FemMultipatchModeler.SetupModelPart may only run once.")

        next_node_id = 1
        next_element_id = 1
        for index, settings in enumerate(self.patch_settings):
            patch = self.model_part.CreateSubModelPart(f"Patch_{index}")
            self._generate_mesh(
                patch, settings, next_node_id, next_element_id, index)
            KM.VariableUtils().SetFlag(KM.ACTIVE, True, patch.Elements)
            next_node_id += patch.NumberOfNodes()
            next_element_id += patch.NumberOfElements()
            self.patches.append(patch)

        background = self.patches[0]
        # Keep the C++ utility alive so background edge topology can be reused
        # by subsequent activation/extraction operations on this fixed mesh.
        self.background_utility = KM.FemMultipatchUtilities(background)
        if self.interface_mode == "gap_sbm":
            inactive_owners = KM.FemMultipatchGapUtilities.UpdateActivation(
                background, self.patches[1:], self.tolerance)
        else:
            inactive_owners = self.background_utility.UpdateBackgroundActivation(
                self.boxes[1:], self.tolerance)
        edges = self.background_utility.ExtractSurrogateEdges(inactive_owners)

        interfaces = self.model_part.CreateSubModelPart("Interfaces")
        next_geometry_id = 1
        for index, bounds in enumerate(self.boxes[1:], 1):
            # Retain per-hole element membership as well as ACTIVE flags.
            inactive = background.CreateSubModelPart(f"Inactive_Patch_{index}")
            element_ids = [
                element_id for element_id, owner in inactive_owners.items()
                if owner == index]
            inactive.AddElements(element_ids)
            inactive.AddNodes(sorted({node.Id for element in inactive.Elements
                                      for node in element.GetGeometry()}))

            interface = interfaces.CreateSubModelPart(f"Patch_{index}")
            if self.interface_mode == "gap_sbm":
                KM.FemMultipatchGapUtilities.CreateGeometry(
                    background, self.patches[index], interface, self.tolerance)
                continue
            true_boundary = interface.CreateSubModelPart("TrueBoundary")
            surrogate = interface.CreateSubModelPart("SurrogateBoundary")
            xmin, ymin, xmax, ymax = bounds
            corners = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
            nodes = [true_boundary.CreateNewNode(next_node_id + i, x, y, 0.0)
                     for i, (x, y) in enumerate(corners)]
            next_node_id += 4
            for i in range(4):
                self._add_line(true_boundary, nodes[i], nodes[(i + 1) % 4],
                               next_geometry_id)
                next_geometry_id += 1

            # Extracted independently: the true contour is never used here.
            for first_id, second_id in edges[index]:
                first = background.GetNode(first_id)
                second = background.GetNode(second_id)
                surrogate.AddNodes([first_id, second_id])
                self._add_line(surrogate, first, second, next_geometry_id)
                next_geometry_id += 1

    @staticmethod
    def _add_line(model_part, first, second, geometry_id):
        line = KM.Line2D2(first, second)
        line.SetId(geometry_id)
        model_part.AddGeometry(line)

    def _read_box(self, settings):
        settings = settings.Clone()
        settings.ValidateAndAssignDefaults(KM.Parameters('''{
            "origin": [0.0, 0.0],
            "dimensions": [1.0, 1.0],
            "divisions": [1, 1],
            "element_type": "triangle",
            "element_name": "",
            "rotation_angle": 0.0
        }'''))
        if not math.isfinite(settings["rotation_angle"].GetDouble()):
            raise ValueError("rotation_angle must be finite (degrees).")
        if settings["element_type"].GetString() not in (
                "triangle", "quadrilateral"):
            raise ValueError(
                "element_type must be 'triangle' or 'quadrilateral'.")
        for key in ("origin", "dimensions", "divisions"):
            if settings[key].size() != 2:
                raise ValueError(f"{key} must have two components (x, y).")
        origin = list(settings["origin"].GetVector())
        dimensions = list(settings["dimensions"].GetVector())
        divisions = [v.GetInt() for v in settings["divisions"].values()]
        if not all(math.isfinite(value) for value in origin + dimensions):
            raise ValueError("Box coordinates and dimensions must be finite.")
        if any(value <= self.tolerance for value in dimensions):
            raise ValueError("Box dimensions must be greater than tolerance.")
        if any(value <= 0 for value in divisions):
            raise ValueError("Box divisions must be positive integers.")
        if any(size / count <= self.tolerance
               for size, count in zip(dimensions, divisions)):
            raise ValueError("Mesh spacing must be greater than tolerance.")
        return settings

    @staticmethod
    def _bounds(settings):
        x, y = settings["origin"].GetVector()
        width, height = settings["dimensions"].GetVector()
        return (x, y, x + width, y + height)

    def _align_internal_box(self, patch, background):
        # Snap round-off-sized differences to the exact background grid.
        # Reject cut elements; centroid classification would leave overlaps.
        origin = background["origin"].GetVector()
        dimensions = background["dimensions"].GetVector()
        bounds = self._bounds(patch)
        snapped = list(bounds)
        for axis in range(2):
            spacing = dimensions[axis] / background["divisions"][axis].GetInt()
            for side in (axis, axis + 2):
                coordinate = origin[axis] + round(
                    (bounds[side] - origin[axis]) / spacing) * spacing
                if abs(coordinate - bounds[side]) > self.tolerance:
                    raise ValueError(
                        "Internal boxes must align with "
                        "background grid faces.")
                snapped[side] = coordinate
            patch["origin"][axis].SetDouble(snapped[axis])
            patch["dimensions"][axis].SetDouble(
                snapped[axis + 2] - snapped[axis])

    def _check_internal_box(self, bounds):
        xmin, ymin, xmax, ymax = bounds
        bxmin, bymin, bxmax, bymax = self.boxes[0]
        tol = self.tolerance
        if not (bxmin + tol < xmin < xmax < bxmax - tol and
                bymin + tol < ymin < ymax < bymax - tol):
            raise ValueError(
                "Internal boxes must lie strictly inside the background.")
        for other in self.boxes[1:]:
            if not (xmax < other[0] - tol or xmin > other[2] + tol or
                    ymax < other[1] - tol or ymin > other[3] + tol):
                raise ValueError("Internal boxes must not overlap or touch.")

    @staticmethod
    def _rotated_corners(settings):
        xmin, ymin, xmax, ymax = FemMultipatchModeler._bounds(settings)
        cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
        angle = math.radians(settings["rotation_angle"].GetDouble())
        c, s = math.cos(angle), math.sin(angle)
        return [(cx + c * (x - cx) - s * (y - cy),
                 cy + s * (x - cx) + c * (y - cy))
                for x, y in ((xmin, ymin), (xmin, ymax),
                             (xmax, ymax), (xmax, ymin))]

    @staticmethod
    def _generate_mesh(patch, settings, start_node_id, start_element_id,
                       properties_id):
        # StructuredMeshGeneratorProcess expects clockwise quadrilateral nodes.
        # Its first local direction is consequently y, and the second is x.
        corners = FemMultipatchModeler._rotated_corners(settings)
        geometry = KM.Quadrilateral2D4(*[
            KM.Node(i + 1, x, y, 0.0) for i, (x, y) in enumerate(corners)])
        mesh_settings = KM.Parameters('''{
            "element_name": "Element2D3N",
            "number_of_divisions": [1, 1, 1],
            "create_skin_sub_model_part": false,
            "elements_properties_id": 0,
            "start_node_id": 1,
            "start_element_id": 1
        }''')
        nx, ny = [value.GetInt() for value in settings["divisions"].values()]
        mesh_settings["number_of_divisions"].SetVector(KM.Vector([ny, nx, 1]))
        if settings["element_type"].GetString() == "quadrilateral":
            mesh_settings["element_name"].SetString("Element2D4N")
        if settings["element_name"].GetString():
            mesh_settings["element_name"].SetString(settings["element_name"].GetString())
        mesh_settings["start_node_id"].SetInt(start_node_id)
        mesh_settings["start_element_id"].SetInt(start_element_id)
        mesh_settings["elements_properties_id"].SetInt(properties_id)
        KM.StructuredMeshGeneratorProcess(
            geometry, patch, mesh_settings).Execute()
