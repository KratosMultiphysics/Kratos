# Importing the Kratos Library
import KratosMultiphysics as KM
import numpy as np
from math import comb
from pathlib import Path
import vtk


def Factory(settings, model):
    return IgaBezierVtuOutputProcess(model, settings["Parameters"])


class IgaBezierVtuOutputProcess(KM.OutputProcess):
    def __init__(self, model, params):
        KM.OutputProcess.__init__(self)
        default_parameters = KM.Parameters("""{
            "output_file_name"                    : "",
            "brep_surface_ids"                    : [],
            "model_part_name"                     : "",
            "nodal_solution_step_data_variables"  : []
        }""")
        params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[params["model_part_name"].GetString()]
        self.output_file_name = Path(
            params["output_file_name"].GetString()).with_suffix(".vtu")
        self.brep_surface_ids = [
            params["brep_surface_ids"][i].GetInt()
            for i in range(params["brep_surface_ids"].size())
        ]
        self.nodal_variables = [
            KM.KratosGlobals.GetVariable(
                params["nodal_solution_step_data_variables"][i].GetString())
            for i in range(params["nodal_solution_step_data_variables"].size())
        ]
        self._done = False
        self._bernstein_inv_cache = {}

    def IsOutputStep(self):
        return not self._done

    def PrintOutput(self):
        points = vtk.vtkPoints()
        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.SetPoints(points)

        cell_degrees = []
        point_weights = []
        field_arrays = {var.Name(): [] for var in self.nodal_variables}
        field_comps  = {}

        for brep_id in self.brep_surface_ids:
            brep = self.model_part.GetGeometry(brep_id)
            thb = brep.GetLocalRefinedSurface() \
                if hasattr(brep, "GetLocalRefinedSurface") else brep

            cp_weight = thb.WeightsByNodeId() if hasattr(thb, "WeightsByNodeId") else {}

            level_meta = _build_level_meta(thb)

            if hasattr(brep, "GetActiveCells"):
                cells = brep.GetActiveCells()
            else:
                su = list(brep.KnotsU())
                sv = list(brep.KnotsV())
                cells = [(su[i], su[i+1], sv[j], sv[j+1])
                         for j in range(len(sv)-1) for i in range(len(su)-1)]

            for (u0, u1, v0, v1) in cells:
                if abs(u1 - u0) < 1e-12 or abs(v1 - v0) < 1e-12:
                    continue

                level_degrees = _find_cell_level(u0, u1, v0, v1, level_meta)
                if level_degrees is None:
                    raise RuntimeError(
                        f"Could not locate active cell ({u0}, {u1}) x "
                        f"({v0}, {v1}) in any level's knot spans for brep {brep_id}.")
                p_u, p_v = level_degrees

                ### Sample at (p_u+1) x (p_v+1) nodes on the cell ###
                num_bezier_nodes_u = p_u + 1
                num_bezier_nodes_v = p_v + 1
                bezier_nodes_u = np.linspace(u0, u1, num_bezier_nodes_u)
                bezier_nodes_v = np.linspace(v0, v1, num_bezier_nodes_v)
                sampled_positions   = np.zeros((num_bezier_nodes_u, num_bezier_nodes_v, 3))
                sampled_denominator = np.zeros((num_bezier_nodes_u, num_bezier_nodes_v))
                field_vals = {}

                local_coords = KM.Array3()
                for i in range(num_bezier_nodes_u):
                    for j in range(num_bezier_nodes_v):
                        local_coords[0] = float(bezier_nodes_u[i])
                        local_coords[1] = float(bezier_nodes_v[j])
                        local_coords[2] = 0.0
                        global_coords = brep.GlobalCoordinates(local_coords)
                        sampled_positions[i, j] = [
                            global_coords[0], global_coords[1], global_coords[2]]
                        sampled_denominator[i, j] = _denominator_W(
                            brep, local_coords, cp_weight)
                        for var in self.nodal_variables:
                            val = _eval_variable(brep, self.model_part, local_coords, var)
                            if var.Name() not in field_vals:
                                ncomp = 1 if np.isscalar(val) else len(val)
                                field_vals[var.Name()] = np.zeros((num_bezier_nodes_u, num_bezier_nodes_v, ncomp))
                                field_comps[var.Name()] = ncomp
                            if field_comps[var.Name()] == 1:
                                field_vals[var.Name()][i, j, 0] = float(val)
                            else:
                                for k in range(field_comps[var.Name()]):
                                    field_vals[var.Name()][i, j, k] = float(val[k])

                ###Invert Bernstein to Bezier CPs ###
                Minv_u = self._bernstein_inv(p_u)
                Minv_v = self._bernstein_inv(p_v)
                weights_Bezier = Minv_u @ sampled_denominator @ Minv_v.T
                weighted_points_Bezier = np.empty_like(sampled_positions)
                for c in range(3):
                    weighted_points_Bezier[..., c] = \
                        Minv_u @ (sampled_denominator * sampled_positions[..., c]) @ Minv_v.T
                points_Bezier = weighted_points_Bezier / weights_Bezier[..., np.newaxis]

                field_values_Bezier = {}
                for var_name, sampled_field in field_vals.items():
                    ncomp = field_comps[var_name]
                    field_Bezier = np.empty_like(sampled_field)
                    for c in range(ncomp):
                        weighted_field = sampled_denominator * sampled_field[..., c]
                        field_Bezier[..., c] = \
                            (Minv_u @ weighted_field @ Minv_v.T) / weights_Bezier
                    field_values_Bezier[var_name] = field_Bezier

                ### Output in VTK canonical order ###
                ien = []
                for (i, j) in _tensor_index_to_vtk(p_u, p_v):
                    ien.append(points.InsertNextPoint(
                        float(points_Bezier[i, j, 0]),
                        float(points_Bezier[i, j, 1]),
                        float(points_Bezier[i, j, 2])))
                    point_weights.append(float(weights_Bezier[i, j]))
                    for var_name, field_Bezier in field_values_Bezier.items():
                        if field_comps[var_name] == 1:
                            field_arrays[var_name].append(float(field_Bezier[i, j, 0]))
                        else:
                            field_arrays[var_name].append(
                                tuple(float(x) for x in field_Bezier[i, j, :]))
                ugrid.InsertNextCell(vtk.VTK_BEZIER_QUADRILATERAL, len(ien), ien)
                cell_degrees.append((p_u, p_v, 0))

        n_pts = points.GetNumberOfPoints()

        w = vtk.vtkDoubleArray()
        w.SetName("RationalWeights")
        w.SetNumberOfComponents(1)
        w.SetNumberOfTuples(n_pts)
        for i, W_i in enumerate(point_weights):
            w.SetValue(i, W_i)
        ugrid.GetPointData().SetRationalWeights(w)

        deg = vtk.vtkIntArray()
        deg.SetName("HigherOrderDegrees")
        deg.SetNumberOfComponents(3)
        for (pu, pv, pw) in cell_degrees:
            deg.InsertNextTuple3(pu, pv, pw)
        ugrid.GetCellData().SetHigherOrderDegrees(deg)

        for var in self.nodal_variables:
            name = var.Name()
            ncomp = field_comps.get(name, 1)
            arr = vtk.vtkDoubleArray()
            arr.SetName(name)
            arr.SetNumberOfComponents(ncomp)
            arr.SetNumberOfTuples(n_pts)
            for i, val in enumerate(field_arrays[name]):
                if ncomp == 1:
                    arr.SetTuple1(i, val)
                else:
                    for k in range(ncomp):
                        arr.SetComponent(i, k, val[k])
            ugrid.GetPointData().AddArray(arr)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputData(ugrid)
        writer.SetFileName(str(self.output_file_name))
        writer.SetDataModeToAscii()
        writer.Write()
        self._done = True

    def _bernstein_inv(self, p):
        if p in self._bernstein_inv_cache:
            return self._bernstein_inv_cache[p]
        t = np.linspace(0.0, 1.0, p + 1)
        M = np.zeros((p + 1, p + 1))
        for k in range(p + 1):
            for i in range(p + 1):
                M[k, i] = comb(p, i) * (t[k] ** i) * ((1.0 - t[k]) ** (p - i))
        Minv = np.linalg.inv(M)
        self._bernstein_inv_cache[p] = Minv
        return Minv


def _unique_breakpoints(knots, tol=1e-10):
    out = []
    for k in knots:
        k = float(k)
        if not out or k - out[-1] > tol:
            out.append(k)
    return out


def _build_level_meta(thb):
    meta = []
    for l in range(thb.NumberOfLevels()):
        level = thb.Levels()[l]
        meta.append({
            "p_u": int(level.DegreeU),
            "p_v": int(level.DegreeV),
            "unique_kv_u": _unique_breakpoints(level.KnotsU),
            "unique_kv_v": _unique_breakpoints(level.KnotsV),
        })
    return meta


def _find_cell_level(u0, u1, v0, v1, level_meta, tol=1e-10):
    for l in range(len(level_meta) - 1, -1, -1):
        m = level_meta[l]
        if _has_span(m["unique_kv_u"], u0, u1, tol) and _has_span(m["unique_kv_v"], v0, v1, tol):
            return m["p_u"], m["p_v"]
    return None


def _has_span(unique_kv, a, b, tol):
    for i in range(len(unique_kv) - 1):
        if abs(unique_kv[i] - a) <= tol and abs(unique_kv[i + 1] - b) <= tol:
            return True
    return False


def _denominator_W(brep, local_coords, cp_weight):
    if not cp_weight:
        return 1.0
    ids, N = brep.EvaluateShapeFunctionsAtLocalCoordinates(local_coords, 0)
    denom = 0.0
    for cp_id, N_i in zip(ids, N):
        w_i = cp_weight.get(int(cp_id), 1.0)
        denom += N_i / w_i
    return 1.0 / denom if denom > 0.0 else 1.0


def _eval_variable(brep, mp, local_coords, variable):
    ids, N = brep.EvaluateShapeFunctionsAtLocalCoordinates(local_coords, 0)
    value = mp.GetNode(ids[0]).GetSolutionStepValue(variable)
    if hasattr(value, "__len__"):
        result = np.zeros(len(value))
        for i, nid in enumerate(ids):
            result += N[i] * np.array(mp.GetNode(nid).GetSolutionStepValue(variable))
        return result
    else:
        result = 0.0
        for i, id in enumerate(ids):
            node = mp.model_part.GetNode(id)
            result += N[i] * node.GetSolutionStepValue(variable)
    return result


def _tensor_index_to_vtk(p_u, p_v):
    order = [(0, 0), (p_u, 0), (p_u, p_v), (0, p_v)]
    if p_u < 2 and p_v < 2:
        return order
    for i in range(1, p_u): order.append((i, 0))              # edge 0: +u
    for j in range(1, p_v): order.append((p_u, j))            # edge 1: +v
    for i in range(p_u - 1, 0, -1): order.append((i, p_v))    # edge 2: -u
    for j in range(p_v - 1, 0, -1): order.append((0, j))      # edge 3: -v
    for j in range(1, p_v):                                   # face interior
        for i in range(1, p_u):
            order.append((i, j))
    return order
