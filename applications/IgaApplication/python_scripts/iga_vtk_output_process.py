# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
import numpy as np
from pathlib import Path

VTK_QUAD = 9
VTK_TRIANGLE = 5
VTK_LINE = 3
VERSION = (2, 4)
GLOBAL_TYPE = 'UnstructuredGrid'

try:
    import h5py as h5
except ImportError as e:
    raise ImportError("pip install h5py") from e

def Factory(settings, model):
    return IgaVTKOutputProcess(model, settings["Parameters"])

class IgaVTKOutputProcess(KM.OutputProcess):

    def __init__(self, model, params):
        KM.OutputProcess.__init__(self)

        # Default parameters for the process
        default_parameters = KM.Parameters("""{
            "output_file_name"     : "",
            "brep_surface_ids"     : [],
            "brep_curve_ids": [],
            "model_part_name"      : "",
            "nodal_solution_step_data_variables" : [],
            "output_refinement_surface"    : [],
            "output_refinement_curve"    : [],
            "output_control_type"  : "step",
            "output_frequency"     : 1,
            "output_interval"      : 0.0,
            "interval"               : [0.0, "End"]
        }""")

        params.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.model_part = model[params["model_part_name"].GetString()]
        self.output_file_name = Path(params["output_file_name"].GetString())

        self.settings = params

        # Ensure correct file extension
        self.output_file_name = self.output_file_name.with_suffix(".vtkhdf")

        # Remove old VTKHDF file at the beginning of a fresh analysis
        if self.output_file_name.exists():
            self.output_file_name.unlink()

        # IDs of Brep surfaces to be exported
        self.brep_surface_ids = [
            params["brep_surface_ids"][i].GetInt()
            for i in range(params["brep_surface_ids"].size())
        ]

        # IDs of Brep curves to be exported
        self.brep_curve_ids = [
            params["brep_curve_ids"][i].GetInt()
            for i in range(params["brep_curve_ids"].size())
        ]

        # Refinement level per parametric direction
        self.output_refinement_surface = [
            params["output_refinement_surface"][i].GetInt()
            for i in range(params["output_refinement_surface"].size())
        ]

        # Refinement level per parametric direction
        self.output_refinement_curve = [
            params["output_refinement_curve"][i].GetInt()
            for i in range(params["output_refinement_curve"].size())
        ]

        self.printed_step_count = 0

        # Variables to export (nodal solution step data)
        self.nodal_variables = []
        for i in range(params["nodal_solution_step_data_variables"].size()):
            var_name = params["nodal_solution_step_data_variables"][i].GetString()
            self.nodal_variables.append(KM.KratosGlobals.GetVariable(var_name))

        for var in self.nodal_variables:
            if not self.model_part.HasNodalSolutionStepVariable(var):
                raise Exception(
                    f"Variable '{var.Name()}' is not in the model part. "
                    f"Make sure to call AddNodalSolutionStepVariable."
                )

        # Cached parametric sampling (shared across time steps)
        self.cached_uv = None

        KM.Logger.PrintWarning(
            "IgaVTKOutputProcess",
            "VTKHDF requires ParaView 5.11 or newer. Older versions may not support it reliably."
        )

    def ExecuteBeforeSolutionLoop(self):
        self.__controller = KM.OutputController(self.model, self.settings)

    def IsOutputStep(self):
        return self.__controller.Evaluate()

    # Print the output in the HDFVTK file
    def PrintOutput(self):
        time = float(self.model_part.ProcessInfo[KM.TIME])

        with h5.File(self.output_file_name, "a") as file:
            # Create or retrieve root group
            if "VTKHDF" not in file:
                root = file.create_group("VTKHDF")
                root.attrs["Version"] = VERSION
                root.attrs["Type"] = GLOBAL_TYPE
            else:
                root = file["VTKHDF"]

            # Rebuild parametric sampling if needed when reopening the file
            if self.cached_uv is None and "Points" in root:
                all_local_coordinates = []
                self.cached_uv_sizes = []

                # BRep surfaces
                for brep_id in self.brep_surface_ids:
                    brep_surface = self.model_part.GetGeometry(brep_id)

                    _, _, _, _, uv = self.__compute_full_grid_brep_surface(
                        brep_surface
                    )

                    all_local_coordinates.extend(uv)
                    self.cached_uv_sizes.append(len(uv))

                # IGA beam curves
                for curve_id in self.brep_curve_ids:
                    curve = self.model_part.GetGeometry(curve_id)

                    _, _, _, _, u_coordinates = self.__compute_full_grid_brep_curve(
                        curve
                    )

                    # Store curve coordinates in the same format as surface coordinates
                    all_local_coordinates.extend(
                        [(u, 0.0) for u in u_coordinates]
                    )
                    self.cached_uv_sizes.append(len(u_coordinates))

                self.cached_uv = all_local_coordinates

            # Geometry is written only once
            if "Points" not in root:
                self.cached_uv_sizes = []

                all_pts = []
                all_conn = []
                all_offsets = []
                all_types = []
                all_uv = []

                point_shift = 0
                conn_shift = 0

                # BRep surfaces
                for brep_id in self.brep_surface_ids:
                    brep_surface = self.model_part.GetGeometry(brep_id)

                    pts, conn, offs, types, uv = (
                        self.__compute_full_grid_brep_surface(brep_surface)
                    )

                    if len(pts) == 0:
                        raise RuntimeError(
                            f"No visualization points were generated for "
                            f"BRep surface {brep_id}."
                        )

                    self.cached_uv_sizes.append(len(uv))

                    # Shift indices for global assembly
                    conn_shifted = conn + point_shift
                    offs_shifted = offs[:-1] + conn_shift

                    all_pts.append(pts)
                    all_conn.append(conn_shifted)
                    all_offsets.append(offs_shifted)
                    all_types.append(types)

                    # Surface local coordinates: (u, v)
                    all_uv.extend(uv)

                    point_shift += len(pts)
                    conn_shift += len(conn)

                # IGA curves / beams
                for curve_id in self.brep_curve_ids:
                    curve = self.model_part.GetGeometry(curve_id)

                    pts, conn, offs, types, u_coordinates = (
                        self.__compute_full_grid_brep_curve(curve)
                    )

                    if len(pts) == 0:
                        raise RuntimeError(
                            f"No visualization points were generated for "
                            f"IGA curve {curve_id}."
                        )

                    self.cached_uv_sizes.append(len(u_coordinates))

                    # Shift indices for global assembly
                    conn_shifted = conn + point_shift
                    offs_shifted = offs[:-1] + conn_shift

                    all_pts.append(pts)
                    all_conn.append(conn_shifted)
                    all_offsets.append(offs_shifted)
                    all_types.append(types)

                    # Store curve coordinates using the same (u, v) format
                    all_uv.extend(
                        [(u, 0.0) for u in u_coordinates]
                    )

                    point_shift += len(pts)
                    conn_shift += len(conn)

                if len(all_pts) == 0:
                    raise RuntimeError(
                        "No geometry was provided for VTK output. "
                        f"brep_surface_ids={self.brep_surface_ids}, "
                        f"brep_curve_ids={self.brep_curve_ids}"
                    )

                # Assemble global arrays
                points = np.vstack(all_pts)
                conn = np.concatenate(all_conn)
                offsets_geom = np.concatenate(all_offsets)
                offsets_geom = np.append(offsets_geom, conn_shift)
                types = np.concatenate(all_types)

                self.cached_uv = all_uv

                # Write geometry
                root.create_dataset("Points", data=points)
                root.create_dataset("Connectivity", data=conn)
                root.create_dataset("Offsets", data=offsets_geom)
                root.create_dataset("Types", data=types)

                root.create_dataset("NumberOfPoints", data=(len(points),))
                root.create_dataset("NumberOfCells", data=(len(types),))
                root.create_dataset("NumberOfConnectivityIds", data=(len(conn),))

                pd = root.create_group("PointData")

                # Initial write of nodal variables
                for var in self.nodal_variables:
                    data = []
                    uv_counter = 0
                    patch_index = 0

                    # Evaluate variables on BRep surfaces
                    for brep_id in self.brep_surface_ids:
                        brep_surface = self.model_part.GetGeometry(brep_id)

                        n_local = self.cached_uv_sizes[patch_index]

                        local_uv = self.cached_uv[
                            uv_counter:uv_counter + n_local
                        ]

                        for u, v in local_uv:
                            val = self.__eval_variable(
                                brep_surface,
                                u,
                                v,
                                var
                            )
                            data.append(val)

                        uv_counter += n_local
                        patch_index += 1

                    # Evaluate variables on IGA curves / beams
                    for curve_id in self.brep_curve_ids:
                        curve = self.model_part.GetGeometry(curve_id)

                        n_local = self.cached_uv_sizes[patch_index]

                        local_uv = self.cached_uv[
                            uv_counter:uv_counter + n_local
                        ]

                        for u, _ in local_uv:
                            val = self.__eval_variable(
                                curve,
                                u,
                                0.0,
                                var
                            )
                            data.append(val)

                        uv_counter += n_local
                        patch_index += 1

                    data = np.array(data)

                    # Ensure 2D array: scalar -> (n, 1)
                    if data.ndim == 1:
                        data = data[:, np.newaxis]

                    pd.create_dataset(
                        var.Name(),
                        data=data,
                        maxshape=(None, data.shape[1])
                    )

                n_pts = len(points)

            else:
                # Append new time-step data
                pd = root["PointData"]
                n_pts = int(root["NumberOfPoints"][0])

                for var in self.nodal_variables:
                    data = []
                    uv_counter = 0
                    patch_index = 0

                    # Evaluate variables on BRep surfaces
                    for brep_id in self.brep_surface_ids:
                        brep_surface = self.model_part.GetGeometry(brep_id)

                        n_local = self.cached_uv_sizes[patch_index]

                        local_uv = self.cached_uv[
                            uv_counter:uv_counter + n_local
                        ]

                        for u, v in local_uv:
                            val = self.__eval_variable(
                                brep_surface,
                                u,
                                v,
                                var
                            )
                            data.append(val)

                        uv_counter += n_local
                        patch_index += 1

                    # Evaluate variables on IGA curves / beams
                    for curve_id in self.brep_curve_ids:
                        curve = self.model_part.GetGeometry(curve_id)

                        n_local = self.cached_uv_sizes[patch_index]

                        local_uv = self.cached_uv[
                            uv_counter:uv_counter + n_local
                        ]

                        for u, _ in local_uv:
                            val = self.__eval_variable(
                                curve,
                                u,
                                0.0,
                                var
                            )
                            data.append(val)

                        uv_counter += n_local
                        patch_index += 1

                    data = np.array(data)

                    if data.ndim == 1:
                        data = data[:, np.newaxis]

                    if data.shape[0] != n_pts:
                        raise RuntimeError(
                            "Non-constant number of points across time steps! "
                            f"Expected {n_pts}, obtained {data.shape[0]}."
                        )

                    ds = pd[var.Name()]

                    old = ds.shape[0]
                    ds.resize((old + n_pts, data.shape[1]))
                    ds[old:] = data

            # Time step bookkeeping
            steps = root.require_group("Steps")

            def create_ds(name):
                if name not in steps:
                    steps.create_dataset(name, (0,), maxshape=(None,), dtype='i8')

            if "Values" not in steps:
                steps.create_dataset("Values", (0,), maxshape=(None,), dtype='f8')

            pdo = steps.require_group("PointDataOffsets")
            for var in self.nodal_variables:
                name = var.Name()
                if name not in pdo:
                    pdo.create_dataset(name, (0,), maxshape=(None,), dtype='i8')

            create_ds("PartOffsets")
            create_ds("NumberOfParts")
            create_ds("PointOffsets")
            create_ds("CellOffsets")
            create_ds("ConnectivityIdOffsets")

            values = steps["Values"]

            offsets_dict = {}
            for var in self.nodal_variables:
                offsets_dict[var.Name()] = pdo[var.Name()]

            part_offsets = steps["PartOffsets"]
            num_parts = steps["NumberOfParts"]
            point_offsets = steps["PointOffsets"]
            cell_offsets = steps["CellOffsets"]
            conn_offsets = steps["ConnectivityIdOffsets"]

            n_steps = values.shape[0]

            # Resize datasets for new step
            values.resize((n_steps + 1,))
            for var in self.nodal_variables:
                offsets_dict[var.Name()].resize((n_steps + 1,))
            part_offsets.resize((n_steps + 1,))
            num_parts.resize((n_steps + 1,))
            point_offsets.resize((n_steps + 1,))
            cell_offsets.resize((n_steps + 1,))
            conn_offsets.resize((n_steps + 1,))

            values[n_steps] = time

            # Offsets mark start of each step block in flattened arrays
            for var in self.nodal_variables:
                offsets = offsets_dict[var.Name()]

                if n_steps == 0:
                    offsets[n_steps] = 0
                else:
                    offsets[n_steps] = offsets[n_steps - 1] + n_pts

            # Geometry is static
            part_offsets[n_steps] = 0
            point_offsets[n_steps] = 0
            cell_offsets[n_steps] = 0
            conn_offsets[n_steps] = 0
            num_parts[n_steps] = 1

            self.printed_step_count += 1
            steps.attrs["NSteps"] = self.printed_step_count

        self.__controller.Update()

    # Compute the visualization grid for brep surfaces
    def __compute_full_grid_brep_surface(self, brep_surface):
        knots_u = brep_surface.KnotsU()
        knots_v = brep_surface.KnotsV()

        u_min, u_max = knots_u[0], knots_u[-1]
        v_min, v_max = knots_v[0], knots_v[-1]

        ku = self.__refine(knots_u, self.output_refinement_surface[0])
        kv = self.__refine(knots_v, self.output_refinement_surface[1])

        pts, conn, offs, types = [], [], [], []
        uv_coords = []

        pid = 0
        c = 0

        # Loop over knot spans and triangulate each parametric cell
        for j in range(len(kv) - 1):
            for i in range(len(ku) - 1):

                u0, u1 = ku[i], ku[i + 1]
                v0, v1 = kv[j], kv[j + 1]

                if abs(u1 - u0) < 1e-12 or abs(v1 - v0) < 1e-12:
                    continue

                trimmed, tris = brep_surface.ComputeSpanTriangulationLocalSpace(u0, u1, v0, v1)

                if not trimmed:
                    coords = [(u0,v0),(u1,v0),(u1,v1),(u0,v1)]
                    ids = []

                    for u,v in coords:
                        u = max(u_min, min(u_max, u))
                        v = max(v_min, min(v_max, v))

                        if not np.isfinite(u) or not np.isfinite(v):
                            continue

                        lc = KM.Array3()
                        lc[0], lc[1], lc[2] = u, v, 0.0

                        try:
                            X = brep_surface.GlobalCoordinates(lc)
                        except Exception as e:
                            raise RuntimeError(
                                f"Failed to compute global coordinates for BRep surface {brep_surface.Id} "
                                f"at local coordinates (u, v)=({u}, {v})."
                            ) from e

                        pts.append([X[0], X[1], X[2]])
                        uv_coords.append((u,v))
                        ids.append(pid)
                        pid += 1

                    conn.extend(ids)
                    offs.append(c)
                    c += 4
                    types.append(VTK_QUAD)

                else:
                    for tri in tris:
                        ids = []
                        for k in range(3):
                            u = tri[k, 0]
                            v = tri[k, 1]

                            u = max(u_min, min(u_max, u))
                            v = max(v_min, min(v_max, v))

                            if not np.isfinite(u) or not np.isfinite(v):
                                continue

                            lc = KM.Array3()
                            lc[0], lc[1], lc[2] = u, v, 0.0

                            try:
                                X = brep_surface.GlobalCoordinates(lc)
                            except Exception as e:
                                raise RuntimeError(
                                    f"Failed to compute global coordinates for BRep surface {brep_surface.Id} "
                                    f"at local coordinates (u, v)=({u}, {v})."
                                ) from e

                            pts.append([X[0], X[1], X[2]])
                            uv_coords.append((u,v))
                            ids.append(pid)
                            pid += 1

                        conn.extend(ids)
                        offs.append(c)
                        c += 3
                        types.append(VTK_TRIANGLE)

        offs.append(c)

        return (
            np.array(pts,float),
            np.array(conn,np.int64),
            np.array(offs,np.int64),
            np.array(types,np.uint8),
            uv_coords
        )

    # Compute the visualization grid for a BRep curve
    def __compute_full_grid_brep_curve(self, brep_curve):
        knots = brep_curve.Knots()

        if len(knots) < 2:
            raise RuntimeError(
                f"BRep curve {brep_curve.Id} has fewer than two knots."
            )

        u_min = knots[0]
        u_max = knots[len(knots) - 1]

        if len(self.output_refinement_curve) == 0:
            refinement = 0
        else:
            refinement = self.output_refinement_curve[0]

        refined_knots = self.__refine(knots, refinement)

        pts = []
        conn = []
        offs = []
        types = []
        u_coordinates = []

        point_id = 0
        connectivity_counter = 0

        # Loop over the refined knot spans
        for i in range(len(refined_knots) - 1):
            u0 = refined_knots[i]
            u1 = refined_knots[i + 1]

            # Skip zero-length knot spans
            if abs(u1 - u0) < 1e-12:
                continue

            local_point_ids = []

            # Create both endpoints of the VTK line
            for u in (u0, u1):
                u = max(u_min, min(u_max, u))

                if not np.isfinite(u):
                    raise RuntimeError(
                        f"Non-finite local coordinate u={u} encountered "
                        f"for BRep curve {brep_curve.Id}."
                    )

                local_coordinates = KM.Array3()
                local_coordinates[0] = u
                local_coordinates[1] = 0.0
                local_coordinates[2] = 0.0

                try:
                    global_coordinates = brep_curve.GlobalCoordinates(
                        local_coordinates
                    )
                except Exception as exception:
                    raise RuntimeError(
                        f"Failed to compute global coordinates for "
                        f"BRep curve {brep_curve.Id} at u={u}."
                    ) from exception

                pts.append([
                    global_coordinates[0],
                    global_coordinates[1],
                    global_coordinates[2],
                ])

                u_coordinates.append(u)
                local_point_ids.append(point_id)
                point_id += 1

            offs.append(connectivity_counter)
            conn.extend(local_point_ids)

            connectivity_counter += 2
            types.append(VTK_LINE)

        offs.append(connectivity_counter)

        if len(pts) == 0:
            raise RuntimeError(
                f"No visualization points were generated for "
                f"BRep curve {brep_curve.Id}. "
                f"Knots: {list(knots)}"
            )

        return (
            np.array(pts, dtype=float),
            np.array(conn, dtype=np.int64),
            np.array(offs, dtype=np.int64),
            np.array(types, dtype=np.uint8),
            u_coordinates,
        )

    # Evaluate the desired variable at a local coordinate position
    def __eval_variable(self, brep, u, v, variable):
        lc = KM.Array3()
        lc[0] = u
        lc[1] = v
        lc[2] = 0.0

        ids, N = brep.EvaluateShapeFunctionsAtLocalCoordinates(lc, 0)

        if len(ids) == 0:
            raise RuntimeError(
                f"No active control points found for geometry {brep.Id} "
                f"at local coordinates ({u}, {v})."
            )

        sample_node = self.model_part.GetNode(ids[0])
        value = sample_node.GetSolutionStepValue(variable)

        if hasattr(value, "__len__"):
            size = len(value)
            result = [0.0] * size

            for i, node_id in enumerate(ids):
                node = self.model_part.GetNode(node_id)
                node_value = node.GetSolutionStepValue(variable)
                weight = N[i]
                for k in range(size):
                    result[k] += weight * node_value[k]
        else:
            result = 0.0

            for i, node_id in enumerate(ids):
                node = self.model_part.GetNode(node_id)
                result += N[i] * node.GetSolutionStepValue(variable)

        return result

    # Uniform subdivision of knot spans
    def __refine(self, knots, n):
        if n == 0:
            return knots

        out = []

        for i in range(len(knots)-1):
            a,b = knots[i], knots[i+1]
            for j in range(n+1):
                t = j/(n+1)
                out.append((1-t)*a + t*b)

        out.append(knots[len(knots) - 1])
        return out
