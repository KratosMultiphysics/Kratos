import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from pathlib import Path
import numpy as np
from typing import Dict, Tuple, List, Optional


# -----------------------------------------------------------------------------
# Factory
# -----------------------------------------------------------------------------
def Factory(Model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    """Kratos factory."""
    if type(parameters) != Kratos.Parameters:
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MarchingCubesInterfaceOutputProcess2(Model, parameters["settings"], optimization_problem)


class MarchingCubesInterfaceOutputProcess2(Kratos.OutputProcess):
    """
    Interface extractor using a watertight Marching Cubes variant for interior
    cells and marching squares caps on boundary faces. Includes robust handling
    for φ≈0 and consistent face pairing (asymptotic decider).
    """

    # ---------------------------------------------------------------------
    # Defaults
    # ---------------------------------------------------------------------
    def GetDefaultParameters(self) -> Kratos.Parameters:
        """Return default Kratos parameters."""
        return Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>_<step>",
                "output_path"                 : "Solid_interface",
                "save_output_files_in_folder" : true,
                "output_interval"             : 1,
                "model_part_name"             : "",

                "flip_sign"                   : true,       // if true, multiplies phi by -1 before processing
                "zero_tol"                    : 1e-12,      // classification stabilizer for phi≈0
                "min_triangle_area"           : 0.0,        // absolute area threshold for sliver filtering

                // Debug for capping
                "debug_caps"                  : false,
                "debug_first_n_caps"          : 1400,

                // Debug for interior
                "debug_interior"              : false,
                "debug_first_n_interior"      : 300,
                "debug_elem_ids"              : [],

                // Optional lookup-table files (txt or json). If not given/found, tables are generated programmatically
                "mc_edge_table_path"          : "",
                "mc_tri_table_path"           : ""
            }
            """
        )

    # ---------------------------------------------------------------------
    # Init
    # ---------------------------------------------------------------------
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # IO
        self.file_name = parameters["file_name"].GetString()
        self.output_file_name_prefix = parameters["file_name"].GetString()
        self.output_path = parameters["output_path"].GetString()
        self.output_path_cloud = "point_cloud"
        self.save_output_files_in_folder = parameters["save_output_files_in_folder"].GetBool()
        self.model_part = model[parameters["model_part_name"].GetString()]
        self.optimization_problem = optimization_problem
        self.output_interval = parameters["output_interval"].GetInt()

        # MC + numerics
        self.flip_sign = parameters["flip_sign"].GetBool()
        self.zero_tol = parameters["zero_tol"].GetDouble()
        self.min_triangle_area = parameters["min_triangle_area"].GetDouble()

        # Debug (caps)
        self.debug_caps = parameters["debug_caps"].GetBool()
        self.debug_first_n_caps = parameters["debug_first_n_caps"].GetInt()
        self._caps_debug_count = 0

        # Debug (interior)
        self.debug_interior = parameters["debug_interior"].GetBool()
        self.debug_first_n_interior = parameters["debug_first_n_interior"].GetInt()
        self.debug_elem_ids = set()
        if parameters["debug_elem_ids"].IsArray():
            for i in range(parameters["debug_elem_ids"].size()):
                self.debug_elem_ids.add(parameters["debug_elem_ids"][i].GetInt())
        self._dbg_interior_count = 0

        # Node registries per PrintOutput() call
        self.created_nodes_by_coord: Dict[Tuple[float, float, float], int] = {}
        self.created_nodes_by_edge: Dict[Tuple[int, int], int] = {}

        # MC lookup tables (edge + tri). If not present, generate.
        edge_path = parameters["mc_edge_table_path"].GetString()
        tri_path = parameters["mc_tri_table_path"].GetString()
        self._MC_EDGE_TABLE, self._MC_TRI_TABLE = self._load_mc_tables(edge_path, tri_path)

    # ---------------------------------------------------------------------
    # Step gating
    # ---------------------------------------------------------------------
    def IsOutputStep(self) -> bool:
        """Return True when it's time to write output for current step."""
        return self.optimization_problem.GetStep() % self.output_interval == 0

    # ---------------------------------------------------------------------
    # PrintOutput entry point
    # ---------------------------------------------------------------------
    def PrintOutput(self) -> None:
        """Entry: construct the interface and write it out."""
        if not self.IsOutputStep():
            return

        # Pull the current control field (phi per node)
        for control in self.optimization_problem.GetListOfControls():
            self.control_field = control.GetPhysicalField().Evaluate()
        if self.flip_sign:
            self.control_field = -1.0 * self.control_field

        # Temporary model part for the interface
        model = Kratos.Model()
        mp = model.CreateModelPart("Interface")
        mp.CreateNewProperties(2)

        # Reset per-call registries
        self.created_nodes_by_coord.clear()
        self.created_nodes_by_edge.clear()
        self._caps_debug_count = 0
        self._dbg_interior_count = 0
        self._face_pairing_by_elem = {}  # (elem_id, face_idx) -> "consecutive" | "across"

        # Build surfaces: interior MC + boundary caps
        self._build_interface_surfaces(mp)

        # Saving
        if self.save_output_files_in_folder:
            self.output_path = Path(self.output_path)
            if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
                kratos_utils.DeleteDirectoryIfExisting(str(self.output_path_cloud))
                self.model_part.ProcessInfo[Kratos.IS_RESTARTED] = True
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
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
        vtu_output.PrintOutput(str(self.output_path / output_file_name))
        Kratos.Logger.PrintInfo(
            "InterfaceOutputProcess",
            f"✅ Interface written to {self.output_path / output_file_name}.vtu"
        )

    # Debug gating ---------------------------------------------------------
    def _dbg_e(self, elem_id: int) -> bool:
        """Should we print interior debug for this element?"""
        if not self.debug_interior:
            return False
        if self.debug_first_n_interior and self._dbg_interior_count >= self.debug_first_n_interior:
            return False
        return (not self.debug_elem_ids) or (elem_id in self.debug_elem_ids)

    # ---------------------------------------------------------------------
    # Main builder: interior MC + boundary caps
    # ---------------------------------------------------------------------
    def _build_interface_surfaces(self, mp_interface: Kratos.ModelPart) -> None:
        """Build MC interior and cap surfaces, write to `mp_interface`."""

        # Vertex ordering (Hex8)
        EDGE_VERTS = np.array(
            [
                [0, 1], [1, 2], [2, 3], [3, 0],
                [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4], [1, 5], [2, 6], [3, 7],
            ],
            dtype=int,
        )

        # Face to vertex cyclic order
        FACE_NODES = [
            (0, 1, 2, 3),  # z-
            (4, 5, 6, 7),  # z+
            (0, 1, 5, 4),  # y-
            (1, 2, 6, 5),  # x+
            (2, 3, 7, 6),  # y+
            (3, 0, 4, 7),  # x-
        ]

        # Face to global edge indices (cyclic, matching FACE_NODES)
        FACE_EDGES = [
            (0,  1,  2,  3),   # face 0: (0,1,2,3)
            (4,  5,  6,  7),   # face 1: (4,5,6,7)
            (0,  9,  4,  8),   # face 2: (0,1,5,4)  -> e0,e9,e4,e8
            (1, 10,  5,  9),   # face 3: (1,2,6,5)  -> e1,e10,e5,e9
            (2, 11,  6, 10),   # face 4: (2,3,7,6)  -> e2,e11,e6,e10
            (3,  8,  7, 11),   # face 5: (3,0,4,7)  -> e3,e8,e7,e11
        ]

        # Build a quick lookup: is (elem_id, face_idx) a boundary face?
        boundary_face_set = set()
        for _elem_b, _fidx_b, _ in self._collect_boundary_faces_hex8():
            boundary_face_set.add((_elem_b.Id, _fidx_b))

        # One-time assertion: FACE_EDGES cycles match FACE_NODES rings
        if not hasattr(self, "_checked_face_edges"):
            geom = next(iter(self.model_part.Elements)).GetGeometry()
            for f_idx, verts in enumerate(FACE_NODES):
                a, b, c, d = verts
                cyc = FACE_EDGES[f_idx]
                ring = [(a, b), (b, c), (c, d), (d, a)]
                ed = [
                    (0, 1), (1, 2), (2, 3), (3, 0),
                    (4, 5), (5, 6), (6, 7), (7, 4),
                    (0, 4), (1, 5), (2, 6), (3, 7),
                ]
                for k, (u, v) in enumerate(ring):
                    assert (u, v) in ed or (v, u) in ed, f"Bad ring {f_idx}"
                    e = cyc[k]
                    eu, ev = ed[e]
                    assert {u, v} == {eu, ev}, f"FACE_EDGES mismatch on face {f_idx}"
            self._checked_face_edges = True

        # Small utils ------------------------------------------------------
        def tri_area(p0, p1, p2) -> float:
            a = np.asarray(p0, float)
            b = np.asarray(p1, float)
            c = np.asarray(p2, float)
            return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))

        def get_or_create_edge_node(edge_key: Tuple[int, int], p: np.ndarray) -> int:
            """Deduplicate by mesh edge id; also register by rounded coordinate so caps can reuse."""
            if edge_key[0] > edge_key[1]:
                edge_key = (edge_key[1], edge_key[0])
            nid = self.created_nodes_by_edge.get(edge_key)
            if nid is not None:
                return nid
            nid = len(self.created_nodes_by_coord) + 1
            mp_interface.CreateNewNode(nid, float(p[0]), float(p[1]), float(p[2]))
            self.created_nodes_by_edge[edge_key] = nid
            key_p = tuple(np.round(p, 12))
            if key_p not in self.created_nodes_by_coord:
                self.created_nodes_by_coord[key_p] = nid
            return nid

        def get_or_create_point_node(p: np.ndarray) -> int:
            key = tuple(np.round(p, 12))
            nid = self.created_nodes_by_coord.get(key)
            if nid is not None:
                return nid
            nid = len(self.created_nodes_by_coord) + 1
            mp_interface.CreateNewNode(nid, float(p[0]), float(p[1]), float(p[2]))
            self.created_nodes_by_coord[key] = nid
            return nid

        def interp_zero(p0, p1, f0, f1) -> np.ndarray:
            """Robust interpolation of the φ=0 point on segment (p0, p1)."""
            eps = 1e-15
            if abs(f0) < eps and abs(f1) < eps:
                t = 0.5
            elif abs(f0) < eps:
                return np.asarray(p0, float)
            elif abs(f1) < eps:
                return np.asarray(p1, float)
            else:
                denom = (f0 - f1)
                t = 0.5 if abs(denom) < eps else f0 / denom
            return np.asarray(p0, float) + t * (np.asarray(p1, float) - np.asarray(p0, float))

        cond_id = mp_interface.NumberOfConditions()
        prop = mp_interface.GetProperties(2)

        # ---- 1) Interior: element-wise MC using lookup tables ----
        for elem in self.model_part.Elements:
            geom = elem.GetGeometry()

            # element scale for tolerances
            all_xyz = np.array([[n.X, n.Y, n.Z] for n in geom], dtype=float)
            h_elem = np.max(np.linalg.norm(all_xyz - all_xyz.mean(axis=0), axis=1))
            plane_tol = max(1e-12, 1e-8 * h_elem)

            # precompute planes (n, d) for the 6 faces; n is unit normal, d s.t. n·x + d = 0 on the face
            face_planes = []
            for f_idx, local in enumerate(FACE_NODES):
                P = np.array([[geom[i].X, geom[i].Y, geom[i].Z] for i in local], dtype=float)
                n = np.cross(P[1] - P[0], P[2] - P[0])
                nn = np.linalg.norm(n)
                if nn < 1e-30:
                    face_planes.append((None, None))  # degenerate; shouldn't happen
                else:
                    n /= nn
                    d = -float(np.dot(n, P[0]))
                    face_planes.append((n, d))

            # raw nodal field (for interpolation & decider)
            phi_raw = np.array([self.control_field[node.Id - 1] for node in geom], dtype=float)

            # biased copy for classification only (prevents flips on φ≈0)
            eps = float(self.zero_tol) if self.zero_tol > 0.0 else 0.0
            bias_sign = -1.0 if np.sum(phi_raw) < 0.0 else 1.0
            phi_cls = phi_raw.copy()
            phi_cls[np.abs(phi_cls) < eps] = eps * (-1.0 if bias_sign < 0.0 else 1.0)

            cubeindex = 0
            for i in range(8):
                if phi_cls[i] < 0.0:
                    cubeindex |= (1 << i)

            edge_mask = int(self._MC_EDGE_TABLE[cubeindex])
            if edge_mask == 0:
                continue

            if self._dbg_e(elem.Id):
                Kratos.Logger.PrintInfo(
                    "InteriorMC",
                    f"elem={elem.Id} cubeindex={cubeindex} edge_mask=0x{edge_mask:03x} "
                    f"phi_raw={np.round(phi_raw,12).tolist()}"
                )

            # Compute all used edge intersections and create/reuse nodes
            edge_points: List[Optional[np.ndarray]] = [None] * 12
            edge_node_ids: List[Optional[int]] = [None] * 12
            for e in range(12):
                if edge_mask & (1 << e):
                    v0, v1 = EDGE_VERTS[e]
                    p0 = np.array([geom[v0].X, geom[v0].Y, geom[v0].Z], dtype=float)
                    p1 = np.array([geom[v1].X, geom[v1].Y, geom[v1].Z], dtype=float)
                    f0 = float(phi_raw[v0])
                    f1 = float(phi_raw[v1])
                    p = interp_zero(p0, p1, f0, f1)
                    edge_points[e] = p
                    g0 = geom[v0].Id
                    g1 = geom[v1].Id
                    nid = get_or_create_edge_node((min(g0, g1), max(g0, g1)), p)
                    edge_node_ids[e] = nid

            # ---- Emit triangles WITHOUT TRI_TABLE (Lewiner-style, with face-labeled adjacency) ----
            def asymptotic_inside_face(f00, f10, f01, f11, eps_dec=1e-30) -> bool:
                a = f00 - f01 - f10 + f11
                if abs(a) < eps_dec:
                    return (0.25 * (f00 + f10 + f01 + f11) < 0.0)
                u = (f00 - f01) / a
                v = (f00 - f10) / a
                if 0.0 <= u <= 1.0 and 0.0 <= v <= 1.0:
                    val = a * u * v + (f10 - f00) * u + (f01 - f00) * v + f00
                else:
                    val = 0.25 * (f00 + f10 + f01 + f11)
                return (val < 0.0)

            # adjacency: edge -> list of (neighbor_edge, face_id)
            neighbors: Dict[int, List[Tuple[int, int]]] = {e: [] for e in range(12) if (edge_mask & (1 << e))}

            def add_pair(ea: int, eb: int, f_id: int) -> None:
                if (edge_mask & (1 << ea)) and (edge_mask & (1 << eb)):
                    neighbors[ea].append((eb, f_id))
                    neighbors[eb].append((ea, f_id))

            for f_idx, face_edges in enumerate(FACE_EDGES):
                fv = FACE_NODES[f_idx]
                Gf = phi_raw[list(fv)]  # RAW φ for decider

                # Skip pure zero faces to avoid duplicate sheets on that face
                face_eps = max(1e-14, 1e-12 * float(np.max(np.abs(Gf)) or 1.0))
                if np.all(np.abs(Gf) <= face_eps):
                    continue

                fcuts = [e for e in face_edges if (edge_mask & (1 << e))]
                if len(fcuts) == 2:
                    add_pair(fcuts[0], fcuts[1], f_idx)

                elif len(fcuts) == 4:
                    # Asymptotic decider -> choose same pairing as MC
                    # Map to bilinear (f00,f10,f01,f11) using order (0,1,3,2)
                    f00, f10, f01, f11 = Gf[0], Gf[1], Gf[3], Gf[2]
                    center_in = asymptotic_inside_face(f00, f10, f01, f11)
                    e0, e1, e2, e3 = face_edges

                    if self._dbg_e(elem.Id):
                        Kratos.Logger.PrintInfo(
                            "InteriorDecider",
                            f"elem={elem.Id} face={f_idx} four_cuts=True "
                            f"diag={'13' if center_in else '02'} center_in={center_in} "
                            f"phi_face={np.round(Gf,12).tolist()}"
                        )
                    if center_in:
                        pairs = [(e0, e1), (e2, e3)]  # consecutive pairing
                        self._face_pairing_by_elem[(elem.Id, f_idx)] = "consecutive"
                    else:
                        pairs = [(e0, e3), (e1, e2)]  # across-the-face pairing
                        self._face_pairing_by_elem[(elem.Id, f_idx)] = "across"
                    for a, b in pairs:
                        add_pair(a, b, f_idx)

            # Walk loops; at each node choose the segment from the *other* incident face
            visited_any: set[int] = set()
            for start in list(neighbors.keys()):
                if start in visited_any or len(neighbors[start]) == 0:
                    continue

                # pick an initial neighbor and remember its face
                nb0, f0 = neighbors[start][0]
                loop_edges = [start, nb0]
                loop_faces = [f0]  # face used to go from [start]->nb0
                prev = start
                cur = nb0
                last_face = f0

                # traverse
                while True:
                    options = neighbors.get(cur, [])
                    nxt = None
                    next_face = None
                    for nb, fid in options:
                        # prefer the neighbor that uses a different face than the one we came from
                        if nb != prev and fid != last_face:
                            nxt = nb
                            next_face = fid
                            break
                    if nxt is None:
                        # fallback: just take the other neighbor (if exists)
                        for nb, fid in options:
                            if nb != prev:
                                nxt = nb
                                next_face = fid
                                break
                    if nxt is None:
                        break  # open path (shouldn't happen in MC)

                    loop_edges.append(nxt)
                    loop_faces.append(next_face)
                    if nxt == start:
                        break
                    prev, cur = cur, nxt
                    last_face = next_face

                # finalize polygon edges (drop repeated start if present)
                edges_in_loop = loop_edges[:-1] if loop_edges and loop_edges[-1] == start else loop_edges
                faces_in_loop = loop_faces[: len(edges_in_loop) - 1] if loop_faces else []

                if len(edges_in_loop) < 3:
                    continue

                # mark ALL edges in this loop as visited to avoid duplicate polygons
                for ee in edges_in_loop:
                    visited_any.add(ee)

                # build vertex arrays from edges, collapsing accidental duplicates
                pts: List[np.ndarray] = []
                ids: List[int] = []
                for e in edges_in_loop:
                    p = edge_points[e]
                    n = edge_node_ids[e]
                    if p is None or n is None:
                        continue
                    if not pts or not np.allclose(p, pts[-1]):
                        pts.append(p)
                        ids.append(n)

                if len(pts) < 3:
                    continue

                # DEBUG: dump the polygon before triangulating
                if self._dbg_e(elem.Id):
                    Kratos.Logger.PrintInfo(
                        "InteriorLoop",
                        f"elem={elem.Id} edges={edges_in_loop} faces={faces_in_loop} "
                        f"nverts={len(pts)}"
                    )
                    for vi, (pp, nn) in enumerate(zip(pts, ids)):
                        Kratos.Logger.PrintInfo(
                            "InteriorLoop",
                            f"  v{vi}: nid={nn} xyz={np.round(pp,12).tolist()}"
                        )

                # robust triangulation of a (possibly non-planar) polygon
                tri_idx = self._triangulate_polygon_3d(pts)  # -> list[(i0,i1,i2)]
                for (i0, i1, i2) in tri_idx:
                    p0, p1, p2 = pts[i0], pts[i1], pts[i2]
                    tri_ids = [ids[i0], ids[i1], ids[i2]]

                    # (a) area filter (optional)
                    if self.min_triangle_area > 0.0:
                        a2 = 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))
                    else:
                        a2 = 1.0
                    if self.min_triangle_area > 0.0 and a2 <= self.min_triangle_area:
                        continue

                    # (b) cull triangles that lie on a *boundary* face plane of this element
                    cull = False
                    for f_idx, (n, d) in enumerate(face_planes):
                        if n is None:
                            continue
                        if (elem.Id, f_idx) not in boundary_face_set:
                            continue
                        r0 = abs(float(np.dot(n, p0) + d))
                        r1 = abs(float(np.dot(n, p1) + d))
                        r2 = abs(float(np.dot(n, p2) + d))
                        if max(r0, r1, r2) <= plane_tol:
                            # interior triangle is coplanar with a boundary face -> let the cap own it
                            if self.debug_caps:
                                Kratos.Logger.PrintInfo(
                                    "InteriorCull",
                                    f"elem={elem.Id} face={f_idx} culled_tri_ids={tri_ids} "
                                    f"r=[{r0:.2e},{r1:.2e},{r2:.2e}] tol={plane_tol:.2e}"
                                )
                            cull = True
                            break
                    if cull:
                        continue

                    # (c) keep it
                    mp_interface.CreateNewCondition("SurfaceCondition3D3N", cond_id, tri_ids, prop)
                    cond_id += 1

        # ---- 2) Boundary caps: marching-squares per exposed face ----
        for elem, face_idx, face_node_ids in self._collect_boundary_faces_hex8():
            geom = elem.GetGeometry()

            # 4 local vertex indices in cyclic order for this face
            local_face = FACE_NODES[face_idx]

            # Coordinates and phi (signed)
            P4 = np.array([[geom[i].X, geom[i].Y, geom[i].Z] for i in local_face], dtype=float)
            F4 = np.array([self.control_field[geom[i].Id - 1] for i in local_face], dtype=float)
            F8 = np.array([self.control_field[n.Id - 1] for n in geom], dtype=float)  # for bias + edge φ

            # element size for snapping tol
            all_xyz = np.array([[n.X, n.Y, n.Z] for n in geom], dtype=float)
            h = np.max(np.linalg.norm(all_xyz - all_xyz.mean(axis=0), axis=1))
            snap_tol = max(1e-12, 1e-8 * h)

            # Element-wide bias sign to resolve near-zero faces consistently (info only)
            bias_sign = -1.0 if np.sum(F8) < 0.0 else 1.0

            # build cap + get debug payload
            pairing = self._face_pairing_by_elem.get((elem.Id, face_idx))  # "consecutive" | "across" | None
            cap_tris, dbg = self._cap_face_from_edges(P4, F4, bias_sign, return_debug=True, pairing_override=pairing)

            if not cap_tris:
                continue

            if self.debug_caps and (self.debug_first_n_caps == 0 or self._caps_debug_count < self.debug_first_n_caps):
                face_globals = [geom[i].Id for i in local_face]
                Kratos.Logger.PrintInfo(
                    "CapFace",
                    f"elem={elem.Id} face_idx={face_idx} local_face={local_face} "
                    f"node_ids={face_globals} phi={np.round(F4,12).tolist()} "
                    f"case={dbg['case']} inside={dbg['inside']} "
                    f"cut_flags={dbg['cut_flags']} cut_pts={dbg['cut_points']} "
                    f"n_tris={dbg['n_tris']} diag={dbg.get('diag','')} "
                    f"override_used={dbg.get('override_used', False)}"
                )

            face_edge_ids = FACE_EDGES[face_idx]

            def snap_point_to_edge_node(pt: np.ndarray) -> Optional[int]:
                """
                Reuse an existing interior edge-node only if 'pt' actually lies
                on that mesh edge (within snap_tol). Otherwise return None.
                """
                # 1) exact coord reuse if created earlier
                key = tuple(np.round(pt, 12))
                nid = self.created_nodes_by_coord.get(key)
                if nid is not None:
                    return nid

                # 2) check each edge of this face: if we have an edge-node for it,
                #    compute the expected interpolation point and compare to 'pt'
                for e_idx in face_edge_ids:
                    v0, v1 = EDGE_VERTS[e_idx]
                    g0 = geom[v0].Id
                    g1 = geom[v1].Id
                    key_edge = (min(g0, g1), max(g0, g1))
                    nid_edge = self.created_nodes_by_edge.get(key_edge)
                    if nid_edge is None:
                        continue
                    # expected intersection on this edge using element φ
                    p0 = np.array([geom[v0].X, geom[v0].Y, geom[v0].Z], dtype=float)
                    p1 = np.array([geom[v1].X, geom[v1].Y, geom[v1].Z], dtype=float)
                    f0 = float(F8[v0])
                    f1 = float(F8[v1])
                    if abs(f0) < 1e-15 and abs(f1) < 1e-15:
                        p_expected = 0.5 * (p0 + p1)
                    elif abs(f0) < 1e-15:
                        p_expected = p0
                    elif abs(f1) < 1e-15:
                        p_expected = p1
                    else:
                        p_expected = p0 + (f0 / (f0 - f1)) * (p1 - p0)
                    if np.linalg.norm(pt - p_expected) <= snap_tol:
                        return nid_edge
                return None

            for t_idx, tri in enumerate(cap_tris):
                p0, p1, p2 = tri[0], tri[1], tri[2]
                if self.min_triangle_area > 0.0 and tri_area(p0, p1, p2) <= self.min_triangle_area:
                    continue
                ids = []
                for p in (p0, p1, p2):
                    nid = snap_point_to_edge_node(p)
                    if nid is None:
                        nid = get_or_create_point_node(p)
                    ids.append(nid)
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", cond_id, ids, prop)

                if self.debug_caps and (self.debug_first_n_caps == 0 or self._caps_debug_count < self.debug_first_n_caps):
                    Kratos.Logger.PrintInfo(
                        "CapFace",
                        f"  tri{t_idx}: node_ids={ids} "
                        f"coords={np.round(p0,12).tolist()} | {np.round(p1,12).tolist()} | {np.round(p2,12).tolist()}"
                    )
                cond_id += 1

            if self.debug_caps and (self.debug_first_n_caps == 0 or self._caps_debug_count < self.debug_first_n_caps):
                self._caps_debug_count += 1

    # ---------------------------------------------------------------------
    # Boundary face collector (unique faces only)
    # ---------------------------------------------------------------------
    def _collect_boundary_faces_hex8(self):
        """
        Yield (elem, face_index, face_node_ids) for faces that appear only once
        in the mesh.
        """
        FACE_NODES = [
            (0, 1, 2, 3),  # z-
            (4, 5, 6, 7),  # z+
            (0, 1, 5, 4),  # y-
            (1, 2, 6, 5),  # x+
            (2, 3, 7, 6),  # y+
            (3, 0, 4, 7),  # x-
        ]
        # map: sorted tuple of 4 global node ids -> (elem, face_idx, local_face_ids)
        face_map: Dict[Tuple[int, int, int, int], Tuple[Kratos.Element, int, Tuple[int, int, int, int]]] = {}
        for elem in self.model_part.Elements:
            geom = elem.GetGeometry()
            for face_idx, local in enumerate(FACE_NODES):
                ids = tuple(sorted([geom[i].Id for i in local]))
                if ids not in face_map:
                    face_map[ids] = (elem, face_idx, local)
                else:
                    # seen twice -> interior face; remove to save memory
                    face_map.pop(ids, None)

        # remaining are boundary
        for (_, (elem, face_idx, local)) in face_map.items():
            geom = elem.GetGeometry()
            yield elem, face_idx, [geom[i].Id for i in local]

    # ---------------------------------------------------------------------
    # Cap construction: marching-squares w/ decider + degeneracy fixes
    # ---------------------------------------------------------------------
    def _cap_face_from_edges(
        self,
        P4: np.ndarray,
        F4s: np.ndarray,
        bias_sign: float = 1.0,
        return_debug: bool = False,
        pairing_override: Optional[str] = None,
    ):
        """
        Marching-squares cap on a quad face (0-1-2-3 cyclic).

        - Uses only edge cut points + existing face corners (no interior points).
        - Asymptotic decider for 5/10 ambiguities; honors `pairing_override`
          ("consecutive" | "across") to match the interior pairing.
        - When center is inside ("consecutive") we keep the 4-tri diamond (watertight).
        - Vertices are snapped to the computed cut points to avoid FP drift.
        """
        P4 = np.asarray(P4, float)
        G = np.asarray(F4s, float).astype(float)

        eps = max(1e-14, 1e-12 * float(np.max(np.abs(G)) or 1.0))

        # all-zero face → fill with 2 tris
        if np.all(np.abs(G) <= eps):
            tris = [
                np.array([P4[0], P4[1], P4[2]], float),
                np.array([P4[0], P4[2], P4[3]], float),
            ]
            if return_debug:
                dbg = {
                    "case": 15,
                    "inside": [True, True, True, True],
                    "phi": np.round(G, 12).tolist(),
                    "cut_flags": [False] * 4,
                    "cut_points": [None] * 4,
                    "n_tris": 2,
                    "diag": "02",
                    "override_used": False,
                }
                return tris, dbg
            return tris

        # classification (same bias as interior) for case id
        Gcls = G.copy()
        if bias_sign == 0.0:
            bias_sign = 1.0
        Gcls[np.abs(Gcls) < eps] = eps * (-1.0 if bias_sign < 0.0 else 1.0)
        inside = (Gcls < 0.0) | np.isclose(Gcls, 0.0)
        b0, b1, b2, b3 = [int(x) for x in inside]
        case = (b0) | (b1 << 1) | (b2 << 2) | (b3 << 3)

        # compute edge cut points with zero-preserving values (exactly on edges)
        pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
        Gz = G.copy()
        Gz[np.abs(Gz) < eps] = 0.0
        cut = [None, None, None, None]
        cut_flags = [False, False, False, False]

        def edge_point(i, j, gi, gj):
            if gi == 0.0 and gj == 0.0:
                return 0.5 * (P4[i] + P4[j])
            if gi == 0.0 and gj != 0.0:
                return P4[i]
            if gj == 0.0 and gi != 0.0:
                return P4[j]
            denom = gi - gj
            t = 0.5 if abs(denom) < 1e-30 else gi / denom
            return P4[i] + t * (P4[j] - P4[i])

        for e, (i, j) in enumerate(pairs):
            gi, gj = float(Gz[i]), float(Gz[j])
            if (gi < 0.0) != (gj < 0.0) or gi == 0.0 or gj == 0.0:
                cut_flags[e] = True
                cut[e] = edge_point(i, j, gi, gj)

        # asymptotic decider on RAW values
        def asymptotic_inside(f00, f10, f01, f11) -> bool:
            a = f00 - f01 - f10 + f11
            if abs(a) < 1e-30:
                return (0.25 * (f00 + f10 + f01 + f11) < 0.0)
            u = (f00 - f01) / a
            v = (f00 - f10) / a
            if 0.0 <= u <= 1.0 and 0.0 <= v <= 1.0:
                val = a * u * v + (f10 - f00) * u + (f01 - f00) * v + f00
            else:
                val = 0.25 * (f00 + f10 + f01 + f11)
            return (val < 0.0)

        f00, f10, f01, f11 = float(G[0]), float(G[1]), float(G[3]), float(G[2])
        center_inside = asymptotic_inside(f00, f10, f01, f11)

        override_used = False
        if case in (5, 10):
            if pairing_override in ("consecutive", "across"):
                use_center = (pairing_override == "consecutive")
                override_used = True
            else:
                use_center = center_inside
        else:
            if pairing_override == "consecutive":
                use_center = True
                override_used = True
            elif pairing_override == "across":
                use_center = False
                override_used = True
            else:
                use_center = False  # default for other cases

        tris: List[np.ndarray] = []

        def add(a, b, c):
            tri = np.array([a, b, c], float)
            if self.min_triangle_area > 0.0:
                area = 0.5 * float(np.linalg.norm(np.cross(tri[1] - tri[0], tri[2] - tri[0])))
                if area <= self.min_triangle_area:
                    return
            tris.append(tri)

        # cases — only use corners and the precomputed edge cuts
        if case == 0:
            pass
        elif case == 15:
            add(P4[0], P4[1], P4[2])
            add(P4[0], P4[2], P4[3])

        elif case == 1:
            add(P4[0], cut[3], cut[0])
        elif case == 2:
            add(P4[1], cut[0], cut[1])
        elif case == 4:
            add(P4[2], cut[1], cut[2])
        elif case == 8:
            add(P4[3], cut[2], cut[3])

        elif case == 3:
            add(P4[0], P4[1], cut[1])
            add(P4[0], cut[1], cut[3])
        elif case == 6:
            add(P4[1], P4[2], cut[2])
            add(P4[1], cut[2], cut[0])
        elif case == 12:
            add(P4[2], P4[3], cut[3])
            add(P4[2], cut[3], cut[1])
        elif case == 9:
            add(P4[3], P4[0], cut[0])
            add(P4[3], cut[0], cut[2])

        elif case == 7:
            A, B = cut[2], cut[3]
            if A is None or B is None or np.allclose(A, B):
                add(P4[0], P4[1], P4[2])
            else:
                add(A, P4[2], P4[1])
                add(A, P4[1], P4[0])
                add(A, P4[0], B)
        elif case == 11:
            A, B = cut[1], cut[2]
            if A is None or B is None or np.allclose(A, B):
                add(P4[0], P4[1], P4[3])
            else:
                add(A, P4[1], P4[0])
                add(A, P4[0], P4[3])
                add(A, P4[3], B)
        elif case == 13:
            A, B = cut[0], cut[1]
            if A is None or B is None or np.allclose(A, B):
                add(P4[0], P4[3], P4[2])
            else:
                add(A, P4[0], P4[3])
                add(A, P4[3], P4[2])
                add(A, P4[2], B)
        elif case == 14:
            A, B = cut[3], cut[0]
            if A is None or B is None or np.allclose(A, B):
                add(P4[1], P4[2], P4[3])
            else:
                add(A, P4[3], P4[2])
                add(A, P4[2], P4[1])
                add(A, P4[1], B)

        elif case == 5:
            if use_center:  # consecutive → “diamond”, 4 tris (watertight)
                add(cut[3], P4[0],  cut[0])
                add(cut[3], cut[0], cut[1])
                add(cut[3], cut[1], P4[2])
                add(cut[3], P4[2],  cut[2])
            else:           # across → 2 tris
                add(P4[0], cut[3], cut[0])
                add(P4[2], cut[1], cut[2])

        elif case == 10:
            if use_center:  # consecutive
                add(cut[0], P4[1],  cut[1])
                add(cut[0], cut[1], cut[2])
                add(cut[0], cut[2], P4[3])
                add(cut[0], P4[3],  cut[3])
            else:           # across
                add(P4[1], cut[0], cut[1])
                add(P4[3], cut[2], cut[3])
        else:
            if self.debug_caps:
                Kratos.Logger.PrintWarning("CapFace", f"Unhandled case {case}")

        # snap every vertex that lies on an edge back to the exact cut point
        if tris:
            L = float(np.linalg.norm(P4.max(axis=0) - P4.min(axis=0)))  # face size scale
            snap_eps = max(1e-12 * L, 1e-14)
            cut_np = [None if c is None else np.asarray(c, float) for c in cut]

            def snap_to_cuts(pt):
                best = None
                bestd = 1e30
                for c in cut_np:
                    if c is None:
                        continue
                    d = float(np.linalg.norm(pt - c))
                    if d < bestd:
                        bestd = d
                        best = c
                return best if best is not None and bestd <= snap_eps else pt

            for t in tris:
                t[0] = snap_to_cuts(t[0])
                t[1] = snap_to_cuts(t[1])
                t[2] = snap_to_cuts(t[2])

        if not return_debug:
            return tris

        dbg = {
            "case": case,
            "inside": inside.tolist(),
            "phi": np.round(G, 12).tolist(),
            "cut_flags": cut_flags,
            "cut_points": [None if c is None else np.round(c, 12).tolist() for c in cut],
            "n_tris": len(tris),
            "diag": ("13" if (case in (5, 10) and use_center) else "02"),
            "override_used": override_used,
        }
        return tris, dbg

    # ---------------------------------------------------------------------
    # Lookup-table loading (txt/json) with programmatic fallback
    # ---------------------------------------------------------------------
    def _load_mc_tables(self, edge_path: str, tri_path: str):
        """Load or generate EDGE_TABLE and TRI_TABLE for MC."""
        here = Path(__file__).resolve().parent

        # Candidate files in order
        edge_candidates = []
        tri_candidates = []
        if edge_path:
            edge_candidates.append(Path(edge_path))
        if tri_path:
            tri_candidates.append(Path(tri_path))
        # defaults next to this file
        edge_candidates += [
            here / "mc_edge_table_generated.txt",
            here / "mc_edge_table.txt",
            here / "mc_edge_table_generated.json",
            here / "mc_edge_table.json",
        ]
        tri_candidates += [
            here / "mc_tri_table_generated.txt",
            here / "mc_tri_table.txt",
            here / "mc_tri_table_generated.json",
            here / "mc_tri_table.json",
        ]

        edge = None
        tri = None

        def read_txt_ints(p: Path) -> Optional[np.ndarray]:
            try:
                arr = np.loadtxt(str(p), dtype=np.int32)
                return arr
            except Exception:
                return None

        def read_json_ints(p: Path) -> Optional[np.ndarray]:
            try:
                import json
                with open(p, "r") as f:
                    data = json.load(f)
                return np.asarray(data, dtype=np.int32)
            except Exception:
                return None

        # Try to read edge table
        for p in edge_candidates:
            arr = read_json_ints(p) if p.suffix.lower() == ".json" else read_txt_ints(p)
            if arr is not None:
                edge = arr
                break

        # Try to read tri table
        for p in tri_candidates:
            arr = read_json_ints(p) if p.suffix.lower() == ".json" else read_txt_ints(p)
            if arr is not None:
                tri = arr
                break

        # Reshape/validate if read
        if edge is not None:
            edge = np.asarray(edge, dtype=np.int32).reshape(256,)
        if tri is not None:
            tri = np.asarray(tri, dtype=np.int32).reshape(256, 16)

        # Generate programmatically if missing (watertight; consistent pairing)
        if edge is None or tri is None:
            edge, tri = self._generate_mc_tables()
            # Persist next to file for future runs
            try:
                np.savetxt(here / "mc_edge_table_generated.txt", edge, fmt="%d")
                np.savetxt(here / "mc_tri_table_generated.txt", tri, fmt="%d")
            except Exception:
                pass

        # Final sanity
        if edge.shape != (256,) or tri.shape != (256, 16):
            raise RuntimeError(f"Invalid MC tables: edge {edge.shape}, tri {tri.shape}")
        return edge, tri

    # ---------------------------------------------------------------------
    # Programmatic generation of EDGE_TABLE and TRI_TABLE (consistent pairing)
    # ---------------------------------------------------------------------
    def _generate_mc_tables(self):
        """Generate MC lookup tables with consistent face pairing."""
        EDGES = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]
        FACES = [
            [0, 1, 2, 3],     # z-
            [4, 5, 6, 7],     # z+
            [3, 11, 7, 8],    # x-
            [1, 10, 5, 9],    # x+
            [2, 10, 6, 11],   # y+
            [0, 9, 4, 8],     # y-
        ]

        EDGE_TABLE = np.zeros(256, dtype=np.int32)
        TRI_TABLE = -np.ones((256, 16), dtype=np.int32)

        for case in range(256):
            inside = [(case >> i) & 1 == 1 for i in range(8)]
            cut_edges = []
            mask = 0
            for e, (a, b) in enumerate(EDGES):
                if inside[a] != inside[b]:
                    mask |= (1 << e)
                    cut_edges.append(e)
            EDGE_TABLE[case] = mask
            if not cut_edges:
                continue

            # For each face, pair the cut edges. If 4 cuts, pair (0-1) and (2-3) in face order.
            neighbors: Dict[int, List[int]] = {e: [] for e in cut_edges}
            for face in FACES:
                fcuts = [e for e in face if e in neighbors]
                if len(fcuts) == 2:
                    a, b = fcuts
                    neighbors[a].append(b)
                    neighbors[b].append(a)
                elif len(fcuts) == 4:
                    a, b, c, d = fcuts
                    for u, v in [(a, b), (c, d)]:
                        neighbors[u].append(v)
                        neighbors[v].append(u)

            # Extract loops and fan-triangulate
            tris: List[Tuple[int, int, int]] = []
            visited = set()
            for start in cut_edges:
                if start in visited or len(neighbors.get(start, [])) == 0:
                    continue
                loop = [start]
                visited.add(start)
                prev = None
                cur = start
                nbrs = neighbors[cur]
                nxt = nbrs[0]
                while True:
                    loop.append(nxt)
                    visited.add(nxt)
                    prev, cur = cur, nxt
                    nbrs = neighbors.get(cur, [])
                    if not nbrs:
                        break
                    if len(nbrs) == 1:
                        nxt = nbrs[0]
                    else:
                        nxt = nbrs[0] if nbrs[1] == prev else nbrs[1]
                    if nxt == start:
                        break
                if loop[-1] != start and start in neighbors.get(loop[-1], []):
                    loop.append(start)
                verts = loop[:-1] if loop and loop[-1] == start else loop
                # fan from 0
                if len(verts) >= 3:
                    v0 = verts[0]
                    for i in range(1, len(verts) - 1):
                        tris.append((v0, verts[i], verts[i + 1]))

            flat = [x for tri in tris for x in tri]
            flat = flat[:16] + [-1] * (16 - len(flat))
            TRI_TABLE[case, :] = np.array(flat, dtype=np.int32)

        return EDGE_TABLE, TRI_TABLE

    # ---------------------------------------------------------------------
    # Helper: interpolation and clipping
    # ---------------------------------------------------------------------
    def _interp_on_edge(self, p0: np.ndarray, p1: np.ndarray, f0: float, f1: float) -> np.ndarray:
        """Robust linear interpolation of the φ=0 point on segment (p0,p1)."""
        p0 = np.asarray(p0, float)
        p1 = np.asarray(p1, float)
        eps = 1e-30
        if abs(f0) < eps and abs(f1) < eps:
            return 0.5 * (p0 + p1)
        if abs(f0) < eps:
            return p0
        if abs(f1) < eps:
            return p1
        denom = (f0 - f1)
        t = 0.5 if abs(denom) < eps else (f0 / denom)
        # clamp tiny overshoots from FP noise
        if t < -1e-12:
            t = 0.0
        if t > 1.0 + 1e-12:
            t = 1.0
        return p0 + t * (p1 - p0)

    def _clip_triangle_phi_le_zero(self, P3x3: np.ndarray, F3: np.ndarray) -> list:
        """
        Clip a single triangle (3x3 points P3x3 with scalar values F3 at the vertices)
        against the half-space φ <= 0 (using self.zero_tol). Returns list of triangles.
        """
        P = np.asarray(P3x3, float)
        F = np.asarray(F3, float)
        tol = float(self.zero_tol)

        inside = [F[i] <= tol for i in range(3)]
        idx_in = [i for i, s in enumerate(inside) if s]
        idx_out = [i for i, s in enumerate(inside) if not s]

        # all in / all out
        if len(idx_in) == 3:
            return [P.copy()]
        if len(idx_in) == 0:
            return []

        # 1 inside, 2 outside -> one clipped triangle
        if len(idx_in) == 1:
            i = idx_in[0]
            j, k = idx_out
            A = self._interp_on_edge(P[i], P[j], F[i], F[j])
            B = self._interp_on_edge(P[i], P[k], F[i], F[k])
            return [np.array([P[i], A, B], dtype=float)]

        # 2 inside, 1 outside -> quad -> two triangles
        if len(idx_in) == 2:
            i, j = idx_in
            k = idx_out[0]
            A = self._interp_on_edge(P[i], P[k], F[i], F[k])
            B = self._interp_on_edge(P[j], P[k], F[j], F[k])
            # Fan triangulation of polygon [Pi, Pj, B, A]
            return [
                np.array([P[i], P[j], B], dtype=float),
                np.array([P[i], B, A], dtype=float),
            ]

    # ---------------------------------------------------------------------
    # Triangulation of a 3D polygon
    # ---------------------------------------------------------------------
    def _triangulate_polygon_3d(self, pts3: List[np.ndarray]) -> List[Tuple[int, int, int]]:
        """
        Triangulate a simple polygon given as an ordered list of 3D points.
        1) fit best-fit plane by SVD,
        2) project to (u,v),
        3) ear-clipping triangulation in 2D,
        4) return triangles as index triplets into 'pts3'.
        Falls back to a fan if polygon is nearly collinear or ear-clipping stalls.
        """
        P = np.asarray(pts3, float)
        n = len(P)
        if n < 3:
            return []

        # best-fit plane (PCA)
        C = P.mean(axis=0)
        Q = P - C
        try:
            _, _, Vt = np.linalg.svd(Q, full_matrices=False)
            u = Vt[0]
            v = Vt[1]
        except Exception:
            # degenerate fallback
            u = np.array([1.0, 0.0, 0.0])
            v = np.array([0.0, 1.0, 0.0])

        # project to 2D
        U = Q @ u
        V = Q @ v
        poly2 = np.column_stack([U, V])

        # signed area (orientation)
        def signed_area(poly):
            s = 0.0
            m = len(poly)
            for i in range(m):
                x1, y1 = poly[i]
                x2, y2 = poly[(i + 1) % m]
                s += x1 * y2 - x2 * y1
            return 0.5 * s

        orient = np.sign(signed_area(poly2))
        if abs(orient) < 1e-18:
            # nearly collinear: fall back to fan
            return [(0, i, i + 1) for i in range(1, n - 1)]

        # ear clipping in 2D
        def is_convex(i0, i1, i2):
            x0, y0 = poly2[i0]
            x1, y1 = poly2[i1]
            x2, y2 = poly2[i2]
            cross = (x2 - x1) * (y0 - y1) - (y2 - y1) * (x0 - x1)
            return (orient * cross) > 0.0

        def point_in_tri(p, a, b, c):
            # barycentric test
            v0 = c - a
            v1 = b - a
            v2 = p - a
            d00 = v0 @ v0
            d01 = v0 @ v1
            d11 = v1 @ v1
            d20 = v2 @ v0
            d21 = v2 @ v1
            denom = d00 * d11 - d01 * d01
            if abs(denom) < 1e-30:
                return False
            v_b = (d11 * d20 - d01 * d21) / denom
            w_b = (d00 * d21 - d01 * d20) / denom
            u_b = 1.0 - v_b - w_b
            eps_b = -1e-14
            return (u_b >= eps_b) and (v_b >= eps_b) and (w_b >= eps_b)

        Vidx = list(range(n))
        tris: List[Tuple[int, int, int]] = []
        guard = 0
        while len(Vidx) >= 3 and guard < 5_000:
            guard += 1
            ear_found = False
            m = len(Vidx)
            for k in range(m):
                i0 = Vidx[(k - 1) % m]
                i1 = Vidx[k]
                i2 = Vidx[(k + 1) % m]

                if not is_convex(i0, i1, i2):
                    continue

                a = poly2[i0]
                b = poly2[i1]
                c = poly2[i2]
                ok = True
                for j in Vidx:
                    if j in (i0, i1, i2):
                        continue
                    if point_in_tri(poly2[j], a, b, c):
                        ok = False
                        break
                if not ok:
                    continue

                tris.append((i0, i1, i2))
                Vidx.pop(k)
                ear_found = True
                break

            if not ear_found:
                # self-intersection or numerical trouble: fallback to fan
                base = Vidx[0]
                for t in range(1, len(Vidx) - 1):
                    tris.append((base, Vidx[t], Vidx[t + 1]))
                break

        return tris


# End of class
