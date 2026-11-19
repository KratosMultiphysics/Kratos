"""Factor of safety of a dry slope with the simplified Bishop method on the mesh of a model part.

The process only reads the mesh and the material properties; it does not solve anything and only
writes the level set of the critical circle to DISTANCE for visualisation. Pore pressures are not
considered (dry slope). It runs in ExecuteFinalizeSolutionStep of the stage it is added to. Note that
the c-phi reduction process replaces the properties of the elements, so add this process to a stage
without it to use the original strength.

Only 2D continuum elements are used: 3- and 6-noded triangles, and 4-, 8- and 9-noded quadrilaterals
(split into two triangles through their corner nodes).
"""

import json
import math

import numpy as np

import KratosMultiphysics as KM

from KratosMultiphysics.GeoMechanicsApplication import bishop_slope_stability as bishop


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return BishopSlopeStabilityProcess(model, settings["Parameters"])


class BishopSlopeStabilityProcess(KM.Process):
    def __init__(self, model, settings):
        super().__init__()
        default_settings = KM.Parameters("""{
            "help": "Factor of safety of a dry slope with the simplified Bishop method",
            "model_part_name": "",
            "number_of_slices": 30,
            "gauss_points_per_slice": 4,
            "search": {
                "center_x_range": [],
                "center_y_range": [],
                "n_centers_x": 10,
                "n_centers_y": 10,
                "radius_range": [],
                "n_radii": 10
            },
            "refine_with_nelder_mead": false,
            "forbidden_boundary_model_parts": [],
            "material_variables": {
                "cohesion": "GEO_COHESION",
                "friction_angle": "GEO_FRICTION_ANGLE",
                "density_solid": "DENSITY_SOLID",
                "porosity": "POROSITY"
            },
            "gravity": 9.81,
            "tolerance": 1e-4,
            "max_iterations": 50,
            "min_m_alpha": 0.2,
            "output_file": "bishop_slope_stability.json",
            "write_distance": true
        }""")
        settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.settings = settings
        self.model = model
        self.model_part = model.GetModelPart(settings["model_part_name"].GetString())
        self.result = None

    # Kratos interface

    def ExecuteFinalizeSolutionStep(self):
        mesh, triangle_element = self._BuildMesh()
        cohesion, friction_angle, unit_weight, material_id = self._ReadMaterials(triangle_element)
        forbidden = self._ForbiddenNodes()
        surface_breaks = bishop.surface_breakpoints(mesh, forbidden if forbidden.any() else None)

        def column_weight(xs, y_base):
            return bishop.column_weights(mesh, unit_weight, xs, y_base)

        def base_properties(x_mid, y_mid):
            points = np.column_stack([x_mid, y_mid])
            index, _ = bishop.locate_points(mesh, points)
            if np.any(index < 0):
                # A midpoint on the outer boundary: take the closest triangle centroid
                centroids = mesh.coordinates[mesh.triangles].mean(axis=1)
                for k in np.flatnonzero(index < 0):
                    index[k] = np.argmin(np.hypot(*(centroids - points[k]).T))
            return cohesion[index], friction_angle[index]

        def extra_breakpoints(arc):
            return np.concatenate([surface_breaks, bishop.base_breakpoints(mesh, arc, material_id)])

        options = dict(number_of_slices=self.settings["number_of_slices"].GetInt(),
                       column_weight=column_weight, base_properties=base_properties,
                       forbidden_nodes=forbidden if forbidden.any() else None,
                       extra_breakpoints=extra_breakpoints,
                       gauss_points_per_slice=self.settings["gauss_points_per_slice"].GetInt(),
                       tolerance=self.settings["tolerance"].GetDouble(),
                       max_iterations=self.settings["max_iterations"].GetInt(),
                       min_m_alpha=self.settings["min_m_alpha"].GetDouble())

        def evaluate(xc, yc, radius):
            if radius <= 0.0:
                return None
            evaluation = bishop.evaluate_circle(mesh, (xc, yc), radius, **options)
            return None if evaluation is None else evaluation[0].factor_of_safety

        centers_x, centers_y, radii = self._SearchGrid(mesh)
        best_factor, best_circle, valid = bishop.grid_search(evaluate, centers_x, centers_y, radii)
        number_of_circles = len(centers_x) * len(centers_y) * len(radii)
        if best_circle is None:
            KM.Logger.PrintWarning("BishopSlopeStabilityProcess", "No admissible slip circle found in the search grid.")
            self.result = {"factor_of_safety": None, "circles_evaluated": number_of_circles, "circles_admissible": 0}
            self._WriteResult()
            return

        refined = False
        if self.settings["refine_with_nelder_mead"].GetBool():
            try:
                refinement = bishop.refine_with_nelder_mead(evaluate, best_circle)
            except ImportError:
                KM.Logger.PrintWarning("BishopSlopeStabilityProcess", "scipy is not available; no refinement.")
                refinement = None
            if refinement is not None and refinement[0] < best_factor:
                best_factor, best_circle = refinement
                refined = True

        result, slices, arc = bishop.evaluate_circle(mesh, best_circle[:2], best_circle[2], **options)
        self.result = self._ResultDictionary(result, slices, arc, number_of_circles, valid, refined)
        self._WriteResult()
        KM.Logger.PrintInfo("BishopSlopeStabilityProcess",
                            f"Bishop factor of safety = {result.factor_of_safety:.4f} for the circle with center "
                            f"({arc.center[0]:.3f}, {arc.center[1]:.3f}) and radius {arc.radius:.3f}")

        if self.settings["write_distance"].GetBool():
            self._WriteDistance(arc)

    # Mesh and data

    def _BuildMesh(self):
        """Triangles (corner nodes) of the continuum elements and the element of each triangle."""
        node_index = {node.Id: i for i, node in enumerate(self.model_part.Nodes)}
        self._node_ids = np.array([node.Id for node in self.model_part.Nodes])
        coordinates = np.array([[node.X0, node.Y0] for node in self.model_part.Nodes])

        triangles, triangle_element = [], []
        for element in self.model_part.Elements:
            geometry = element.GetGeometry()
            if geometry.LocalSpaceDimension() != 2 or not element.IsActive():
                continue
            ids = [node_index[geometry[i].Id] for i in range(geometry.PointsNumber())]
            if len(ids) in (3, 6):
                corner_sets = [ids[:3]]
            elif len(ids) in (4, 8, 9):
                corner_sets = [[ids[0], ids[1], ids[2]], [ids[0], ids[2], ids[3]]]
            else:
                raise RuntimeError(f"BishopSlopeStabilityProcess: element {element.Id} with {len(ids)} nodes "
                                   f"is not supported")
            for corners in corner_sets:
                triangles.append(corners)
                triangle_element.append(element)
        if not triangles:
            raise RuntimeError("BishopSlopeStabilityProcess: no 2D continuum elements in the model part")
        self._coordinates = coordinates
        return bishop.build_triangle_mesh(coordinates, triangles), triangle_element

    def _GetVariable(self, key):
        return KM.KratosGlobals.GetVariable(self.settings["material_variables"][key].GetString())

    @staticmethod
    def _GetPropertyValue(properties, variable, umat_index_variable=None, default=None):
        if properties.Has(variable):
            return properties[variable]
        if umat_index_variable is not None and properties.Has(umat_index_variable) and \
                properties.Has(KM.KratosGlobals.GetVariable("UMAT_PARAMETERS")):
            umat = properties[KM.KratosGlobals.GetVariable("UMAT_PARAMETERS")]
            return umat[properties[umat_index_variable] - 1]
        if default is not None:
            return default
        raise RuntimeError(f"BishopSlopeStabilityProcess: properties {properties.Id} have no {variable.Name()}")

    def _ReadMaterials(self, triangle_element):
        """c', phi' [rad], dry unit weight (1 - n) rho_s g and a material id per triangle."""
        cohesion_variable, friction_variable = self._GetVariable("cohesion"), self._GetVariable("friction_angle")
        density_solid, porosity_variable = self._GetVariable("density_solid"), self._GetVariable("porosity")
        c_index = KM.KratosGlobals.GetVariable("INDEX_OF_UMAT_C_PARAMETER")
        phi_index = KM.KratosGlobals.GetVariable("INDEX_OF_UMAT_PHI_PARAMETER")
        gravity = self.settings["gravity"].GetDouble()

        per_element = {}
        for element in triangle_element:
            if element.Id in per_element:
                continue
            properties = element.Properties
            porosity = self._GetPropertyValue(properties, porosity_variable, default=0.0)
            unit_weight = (1.0 - porosity) * self._GetPropertyValue(properties, density_solid) * gravity
            per_element[element.Id] = (self._GetPropertyValue(properties, cohesion_variable, c_index),
                                       math.radians(self._GetPropertyValue(properties, friction_variable, phi_index)),
                                       unit_weight, properties.Id)

        values = np.array([per_element[element.Id] for element in triangle_element])
        return values[:, 0], values[:, 1], values[:, 2], values[:, 3].astype(int)

    def _ForbiddenNodes(self):
        forbidden = np.zeros(len(self._coordinates), dtype=bool)
        index = {node_id: i for i, node_id in enumerate(self._node_ids)}
        names = self.settings["forbidden_boundary_model_parts"]
        for i in range(names.size()):
            for node in self.model.GetModelPart(names[i].GetString()).Nodes:
                if node.Id in index:
                    forbidden[index[node.Id]] = True
        return forbidden

    def _SearchGrid(self, mesh):
        """Grid of centers and radii; empty ranges default to a box above the model."""
        x_min, y_min = mesh.coordinates.min(axis=0)
        x_max, y_max = mesh.coordinates.max(axis=0)
        width, height = x_max - x_min, y_max - y_min
        search = self.settings["search"]

        def axis(range_key, count_key, default_range):
            values = search[range_key].GetVector() if search[range_key].size() > 0 else default_range
            if len(values) != 2:
                raise RuntimeError(f"BishopSlopeStabilityProcess: search.{range_key} must have two values")
            return np.linspace(values[0], values[1], search[count_key].GetInt())

        return (axis("center_x_range", "n_centers_x", [x_min, x_max]),
                axis("center_y_range", "n_centers_y", [y_max, y_max + height]),
                axis("radius_range", "n_radii", [0.1 * max(width, height), math.hypot(width, 2.0 * height)]))

    def _ResultDictionary(self, result, slices, arc, number_of_circles, valid, refined):
        return {
            "factor_of_safety": result.factor_of_safety,
            "fellenius_factor_of_safety": result.fellenius_factor_of_safety,
            "iterations": result.iterations,
            "center": arc.center.tolist(),
            "radius": arc.radius,
            "entry_point": arc.entry.tolist(),
            "exit_point": arc.exit.tolist(),
            "x_entry": float(arc.entry[0]),
            "x_exit": float(arc.exit[0]),
            "sliding_direction": "negative x" if result.direction > 0 else "positive x",
            "circles_evaluated": number_of_circles,
            "circles_admissible": valid,
            "refined_with_nelder_mead": refined,
            "slices": [{
                "x_left": float(slices.x_left[k]), "x_right": float(slices.x_right[k]),
                "b": float(slices.width[k]), "alpha_deg": math.degrees(slices.alpha[k]),
                "l": float(slices.base_length[k]), "W": float(slices.weight[k]),
                "c": float(slices.cohesion[k]), "phi_deg": math.degrees(slices.friction_angle[k]),
                "m_alpha": float(result.m_alpha[k])
            } for k in range(len(slices.width))]
        }

    def _WriteResult(self):
        with open(self.settings["output_file"].GetString(), "w") as f:
            json.dump(self.result, f, indent=2)

    def _WriteDistance(self, arc):
        historical = self.model_part.HasNodalSolutionStepVariable(KM.DISTANCE)
        for node in self.model_part.Nodes:
            phi = math.hypot(node.X0 - arc.center[0], node.Y0 - arc.center[1]) - arc.radius
            if historical:
                node.SetSolutionStepValue(KM.DISTANCE, phi)
            else:
                node.SetValue(KM.DISTANCE, phi)
