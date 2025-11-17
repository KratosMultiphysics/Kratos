import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from pathlib import Path
import numpy as np
from skimage.measure import marching_cubes

def Factory(Model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if type(parameters) != Kratos.Parameters:
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MarchingCubesInterfaceOutputProcess(Model, parameters["settings"], optimization_problem)


class MarchingCubesInterfaceOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>_<step>",
                "output_path"                 : "Solid_interface",
                "save_output_files_in_folder" : true,
                "output_interval"             : 1,
                "model_part_name"             : "",
                "cap_domain_boundary"         : true,
                "solid_is_negative"           : false,
                "zero_tol"                    : 1e-12,
                "debug_caps"                  : false,
                "debug_first_n_caps"          : 1500,
                "min_triangle_area"           : 1e-12
            }
            """
        )

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.file_name                   = parameters["file_name"].GetString()
        self.output_file_name_prefix     = parameters["file_name"].GetString()
        self.output_path                 = parameters["output_path"].GetString()
        self.output_path_cloud           = "point_cloud"
        self.save_output_files_in_folder = parameters["save_output_files_in_folder"].GetBool()
        self.model_part                  = model[parameters["model_part_name"].GetString()]
        self.optimization_problem        = optimization_problem
        self.output_interval             = parameters["output_interval"].GetInt()
        self.created_coordinates         = {}  # maps (x,y,z) → node_id
        self.cap_domain_boundary         = parameters["cap_domain_boundary"].GetBool()
        self.solid_is_negative           = parameters["solid_is_negative"].GetBool()
        self.zero_tol                    = parameters["zero_tol"].GetDouble()
        self.debug_caps                  = parameters["debug_caps"].GetBool()
        self.debug_first_n_caps          = parameters["debug_first_n_caps"].GetInt()
        self.min_triangle_area           = parameters["min_triangle_area"].GetDouble()

    # -------------------------------------------------------------------------
    def IsOutputStep(self):
        return self.optimization_problem.GetStep() % self.output_interval == 0

    # -------------------------------------------------------------------------
    def PrintOutput(self):
        if not self.IsOutputStep():
            return

        for control in self.optimization_problem.GetListOfControls():
            self.control_field = control.GetPhysicalField().Evaluate()

        model = Kratos.Model()
        mp = model.CreateModelPart("Interface")
        mp.CreateNewProperties(2)
        self.created_coordinates = {}
        self._caps_debug_count = 0

        self.FindIntersectionPointsAndSurfaces_Hexa3D8N(mp)

        # === saving part ===
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
        Kratos.Logger.PrintInfo("InterfaceOutputProcess", f"✅ Interface written to {self.output_path / output_file_name}.vtu")

    # -------------------------------------------------------------------------
    def FindIntersectionPointsAndSurfaces_Hexa3D8N(self, mp_interface: Kratos.ModelPart):
        """For each Hexa3D8N element, extract zero-level iso-surface."""

        # Marching cube vertex order reference (consistent with Kratos Hexa3D8N)
        local_corner_order = np.array([
            [0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1]
        ])

        def get_or_create_node_id(coord):
            """Avoid duplicate interface nodes."""
            key = tuple(np.round(coord, 12))
            if key not in self.created_coordinates:
                new_id = len(self.created_coordinates) + 1
                mp_interface.CreateNewNode(new_id, *coord)
                self.created_coordinates[key] = new_id
            return self.created_coordinates[key]

        cond_id = mp_interface.NumberOfConditions()
        prop = mp_interface.GetProperties(2)

        elem: Kratos.Element
        for elem in self.model_part.Elements:
            geom = elem.GetGeometry()
            # get nodal φ values and coordinates
            phi_vals = np.array([self.control_field[node.Id - 1] for node in geom])
            coords = np.array([[node.X, node.Y, node.Z] for node in geom])

            # map physical to grid cube shape (2x2x2)
            phi_grid = np.zeros((2,2,2))
            for i, (cx, cy, cz) in enumerate(local_corner_order):
                phi_grid[cx, cy, cz] = phi_vals[i]

            try:
                verts, faces, _, _ = marching_cubes(phi_grid, level=0, spacing=(1,1,1))
            except ValueError:
                continue

            if len(verts) == 0:
                continue

            # map vertices from cube to real coordinates via trilinear interpolation
            def trilinear_map(xi, eta, zeta):
                N = np.array([
                    (1 - xi)*(1 - eta)*(1 - zeta),
                    xi*(1 - eta)*(1 - zeta),
                    xi*eta*(1 - zeta),
                    (1 - xi)*eta*(1 - zeta),
                    (1 - xi)*(1 - eta)*zeta,
                    xi*(1 - eta)*zeta,
                    xi*eta*zeta,
                    (1 - xi)*eta*zeta
                ])
                return N @ coords

            # convert [0,1]^3 → [-1,1]^3 for interpolation
            verts_phys = np.array([trilinear_map(v[0], v[1], v[2]) for v in verts])

            for tri in faces:
                p0 = verts_phys[tri[0]]
                p1 = verts_phys[tri[1]]
                p2 = verts_phys[tri[2]]

                # DROP tiny / degenerate triangles
                if self._tri_area(p0, p1, p2) <= self.min_triangle_area:
                    continue

                node_ids = [get_or_create_node_id(p0),
                            get_or_create_node_id(p1),
                            get_or_create_node_id(p2)]
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)
                cond_id += 1
            
        # Add caps on domain boundary faces so the surface closes
        if self.cap_domain_boundary:
            cond_id = self._add_boundary_caps_hex8(mp_interface, get_or_create_node_id, cond_id, prop)

        # Kratos.Logger.PrintInfo("InterfaceOutputProcess", f"✅ Created {len(self.created_coordinates)} interface nodes and {cond_id-1} surface conditions.")

    def _tri_area(self, p0, p1, p2) -> float:
        """Unsigned area of a 3D triangle (0 if degenerate or bad input)."""
        a = np.asarray(p0, dtype=float).reshape(-1)
        b = np.asarray(p1, dtype=float).reshape(-1)
        c = np.asarray(p2, dtype=float).reshape(-1)
        if a.size != 3 or b.size != 3 or c.size != 3:
            return 0.0
        n = np.cross(b - a, c - a)
        return float(np.linalg.norm(n)) * 0.5



    def _interp_on_edge(self, p0, p1, f0, f1):
        # robust linear interpolation for φ = 0 along a segment (p0,p1)
        denom = (f0 - f1)
        t = 0.5 if abs(denom) < 1e-30 else f0 / denom
        return p0 + t * (p1 - p0)

    def _cap_face_from_edges(self, P4, F4s, tol):
        """
        Cap {phi<=0} on a quad face in cyclic corner order (0-1-2-3),
        consistent with marching-cubes, including:
        - strict iso=0 classification,
        - asymptotic decider for 5/10,
        - k=3 degeneracy fix (two cuts collapse at the outside corner).
        """
        P4 = np.asarray(P4, float)
        G  = np.asarray(F4s, float).copy()

        # ---- MC-consistent classification (iso=0) ----
        eps = max(1e-14, 1e-12 * float(np.max(np.abs(G)) or 1.0))
        G[np.abs(G) < eps] = 0.0
        inside = (G < 0.0)
        b0, b1, b2, b3 = [int(x) for x in inside]
        case = (b0) | (b1 << 1) | (b2 << 2) | (b3 << 3)

        # edges e0=(0-1), e1=(1-2), e2=(2-3), e3=(3-0)
        pairs = [(0,1), (1,2), (2,3), (3,0)]
        cut = [None, None, None, None]
        cut_edge = [False, False, False, False]

        def edge_point(i, j, gi, gj):
            denom = gi - gj
            t = 0.5 if abs(denom) < eps else gi / denom
            return P4[i] + t * (P4[j] - P4[i])

        for e, (i, j) in enumerate(pairs):
            gi, gj = G[i], G[j]
            if (gi < 0.0) != (gj < 0.0):
                cut_edge[e] = True
                cut[e] = edge_point(i, j, gi, gj)
            elif gi == 0.0 and gj != 0.0:
                cut_edge[e] = True
                cut[e] = P4[i]
            elif gj == 0.0 and gi != 0.0:
                cut_edge[e] = True
                cut[e] = P4[j]
            # if gi==gj==0: the whole edge is iso; leave it to MC – no cap here

        tris = []

        def tri_area(p0, p1, p2):
            a = np.asarray(p0, float); b = np.asarray(p1, float); c = np.asarray(p2, float)
            return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))

        area_min = getattr(self, "min_triangle_area", 0.0) or 0.0

        def add_tri(a, b, c):
            tri = np.array([a, b, c], dtype=float)
            if tri_area(tri[0], tri[1], tri[2]) > max(area_min, 1e-20):
                tris.append(tri)

        # helper: “cuts coincide?” for degeneracy at a corner
        def same_pt(u, v):
            if u is None or v is None: return False
            L = np.linalg.norm(P4.max(0) - P4.min(0))
            tol_len = max(1e-12 * (L if L > 0 else 1.0), 1e-14)
            return np.linalg.norm(np.asarray(u) - np.asarray(v)) <= tol_len

        # asymptotic decider for ambiguous 5/10
        phi_center = 0.25 * float(G[0] + G[1] + G[2] + G[3])
        center_inside = (phi_center < 0.0)

        # ---- case table ----
        if case == 0:
            pass
        elif case == 15:
            add_tri(P4[0], P4[1], P4[2])
            add_tri(P4[0], P4[2], P4[3])

        elif case == 1:
            add_tri(P4[0], cut[3], cut[0])
        elif case == 2:
            add_tri(P4[1], cut[0], cut[1])
        elif case == 4:
            add_tri(P4[2], cut[1], cut[2])
        elif case == 8:
            add_tri(P4[3], cut[2], cut[3])

        elif case == 3:      # 0,1 in
            add_tri(P4[0], P4[1], cut[1])
            add_tri(P4[0], cut[1], cut[3])
        elif case == 6:      # 1,2 in
            add_tri(P4[1], P4[2], cut[2])
            add_tri(P4[1], cut[2], cut[0])
        elif case == 12:     # 2,3 in
            add_tri(P4[2], P4[3], cut[3])
            add_tri(P4[2], cut[3], cut[1])
        elif case == 9:      # 3,0 in
            add_tri(P4[3], P4[0], cut[0])
            add_tri(P4[3], cut[0], cut[2])

        # ---- k=3 (one outside) 7/11/13/14  — with degeneracy fallback ----
        elif case == 7:      # 0,1,2 in; 3 out -> cuts e2,e3
            A, B = cut[2], cut[3]
            if same_pt(A, B):
                add_tri(P4[0], P4[1], P4[2])   # fallback: 3 inside corners
            else:
                add_tri(A, P4[2], P4[1])
                add_tri(A, P4[1], P4[0])
                add_tri(A, P4[0], B)
        elif case == 11:     # 0,1,3 in; 2 out -> cuts e1,e2
            A, B = cut[1], cut[2]
            if same_pt(A, B):
                add_tri(P4[0], P4[1], P4[3])
            else:
                add_tri(A, P4[1], P4[0])
                add_tri(A, P4[0], P4[3])
                add_tri(A, P4[3], B)
        elif case == 13:     # 0,2,3 in; 1 out -> cuts e0,e1
            A, B = cut[0], cut[1]
            if same_pt(A, B):
                add_tri(P4[0], P4[3], P4[2])
            else:
                add_tri(A, P4[0], P4[3])
                add_tri(A, P4[3], P4[2])
                add_tri(A, P4[2], B)
        elif case == 14:     # 1,2,3 in; 0 out -> cuts e3,e0
            A, B = cut[3], cut[0]
            if same_pt(A, B):
                add_tri(P4[1], P4[2], P4[3])
            else:
                add_tri(A, P4[3], P4[2])
                add_tri(A, P4[2], P4[1])
                add_tri(A, P4[1], B)

        # ---- ambiguous opposite-corners 5/10 with decider ----
        elif case == 5:      # 0,2 in (all edges cut)
            if center_inside:
                add_tri(cut[3], P4[0],  cut[0])
                add_tri(cut[3], cut[0], cut[1])
                add_tri(cut[3], cut[1], P4[2])
                add_tri(cut[3], P4[2],  cut[2])
            else:
                add_tri(P4[0], cut[3], cut[0])
                add_tri(P4[2], cut[1], cut[2])
        elif case == 10:     # 1,3 in (all edges cut)
            if center_inside:
                add_tri(cut[0], P4[1],  cut[1])
                add_tri(cut[0], cut[1], cut[2])
                add_tri(cut[0], cut[2], P4[3])
                add_tri(cut[0], P4[3],  cut[3])
            else:
                add_tri(P4[1], cut[0], cut[1])
                add_tri(P4[3], cut[2], cut[3])

        else:
            if getattr(self, "debug_caps", False):
                Kratos.Logger.PrintWarning("CapFace", f"Unhandled case {case}")

        # ---- debug (optional) ----
        if getattr(self, "debug_caps", False):
            Kratos.Logger.PrintInfo(
                "CapFace",
                f"case={case:02d} center={0.25*float(G.sum()):.3e} inside={inside.tolist()} "
                f"cuts={[bool(x) for x in cut_edge]} n_tris={len(tris)}"
            )
            if getattr(self, "debug_first_n_caps", 0) == 0 or \
            getattr(self, "_caps_debug_count", 0) < self.debug_first_n_caps:
                for k, tri in enumerate(tris):
                    Kratos.Logger.PrintInfo("CapFace",
                        f"  tri{k}: {tri[0].tolist()} | {tri[1].tolist()} | {tri[2].tolist()}")
                try:
                    self._caps_debug_count += 1
                except Exception:
                    pass

        return tris



    def _clip_triangle_phi_le_zero(self, P3x3, F3):
        """
        Clip a single triangle (3D points P3x3, scalar values F3 at vertices)
        against the half-space φ <= 0 (with tol). Returns list of 3D triangles.
        Piecewise-linear inside the face, so the φ=0 boundary is a line segment.
        """
        tol = self.zero_tol
        inside = [f <= tol for f in F3]
        idx_in  = [i for i, s in enumerate(inside) if s]
        idx_out = [i for i, s in enumerate(inside) if not s]

        if len(idx_in) == 3:
            return [P3x3]  # whole tri
        if len(idx_in) == 0:
            return []      # nothing

        # 1 inside, 2 outside -> one clipped triangle
        if len(idx_in) == 1:
            i = idx_in[0]; j, k = idx_out
            A = self._interp_on_edge(P3x3[i], P3x3[j], F3[i], F3[j])
            B = self._interp_on_edge(P3x3[i], P3x3[k], F3[i], F3[k])
            return [np.array([P3x3[i], A, B])]

        # 2 inside, 1 outside -> a quad, triangulate into two
        if len(idx_in) == 2:
            i, j = idx_in; k = idx_out[0]
            A = self._interp_on_edge(P3x3[i], P3x3[k], F3[i], F3[k])
            B = self._interp_on_edge(P3x3[j], P3x3[k], F3[j], F3[k])
            # polygon order [Pi, Pj, B, A] – fan triangulation
            return [
                np.array([P3x3[i], P3x3[j], B]),
                np.array([P3x3[i], B,       A])
            ]

    def _collect_boundary_faces_hex8(self):
        """
        Returns a list of (elem, face_local_nodes) where the face is on the domain boundary.
        Local node ordering for Hexa3D8N faces:
        0: (0,1,2,3)  bottom z=0
        1: (4,5,6,7)  top    z=1
        2: (0,1,5,4)  y=0
        3: (1,2,6,5)  x=1
        4: (2,3,7,6)  y=1
        5: (3,0,4,7)  x=0
        """
        FACE_NODES = [
            (0,1,2,3), (4,5,6,7),
            (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7)
        ]
        face_map = {}
        for elem in self.model_part.Elements:
            geom = elem.GetGeometry()
            elem_node_ids = [node.Id for node in geom]
            for lf, fn in enumerate(FACE_NODES):
                key = tuple(sorted(elem_node_ids[i] for i in fn))
                # first time we see this face -> store; second time -> interior -> mark None
                if key not in face_map:
                    face_map[key] = (elem, fn)
                else:
                    face_map[key] = None
        boundary = []
        for v in face_map.values():
            if v is not None:
                boundary.append(v)
        return boundary

    def _ambiguous_face_mask(self, inside):
        # True when face has opposite corners inside/outside (0101 or 1010)
        return (inside[0] and inside[2] and not inside[1] and not inside[3]) or \
            (inside[1] and inside[3] and not inside[0] and not inside[2])

    def _choose_face_tris(self, F4_signed):
        """
        Pick which diagonal to use for triangulating a Hex8 face so that
        the cap matches marching-cubes on that face.

        F4_signed: φ values at the 4 face corners after applying sgn
                (so 'inside' means F<=0).
        Returns two triangles as index triplets over [0,1,2,3].
        """
        tol = self.zero_tol
        inside = F4_signed <= tol

        # default diagonal: 0-2
        tris = [(0,1,2), (0,2,3)]

        # Ambiguous (opposite corners). Use asymptotic-decider:
        # evaluate φ at face center with bilinear weights (¼ each).
        if self._ambiguous_face_mask(inside):
            phi_c = 0.25 * np.sum(F4_signed)  # signed center value
            # If center is inside, we must cut through the *outside* corners,
            # i.e. use the 1-3 diagonal. Otherwise keep 0-2.
            if phi_c <= tol:
                tris = [(0,1,3), (1,2,3)]
        return tris

    def _add_boundary_caps_hex8(self, mp_interface, get_or_create_node_id, cond_id, prop):
        sgn = 1.0 if self.solid_is_negative else -1.0
        tol = self.zero_tol

        FACE_NODES = [
            (0,1,2,3), (4,5,6,7),           # z- , z+
            (0,1,5,4), (1,2,6,5),           # y- , x+
            (2,3,7,6), (3,0,4,7)            # y+ , x-
        ]

        for elem, face_nodes in self._collect_boundary_faces_hex8():
            geom = elem.GetGeometry()
            nodes = [geom[i] for i in face_nodes]
            P4 = np.array([[n.X, n.Y, n.Z] for n in nodes])
            F4s = sgn * np.array([self.control_field[n.Id - 1] for n in nodes])

            if self.debug_caps and (self.debug_first_n_caps == 0 or self._caps_debug_count < self.debug_first_n_caps):
                inside_dbg = (F4s <= self.zero_tol)
                pairs = [(0,1),(1,2),(2,3),(3,0)]
                cut_flags = []
                cut_pts = []
                for e,(i,j) in enumerate(pairs):
                    fi, fj = F4s[i], F4s[j]
                    chg = (fi <= self.zero_tol) != (fj <= self.zero_tol)
                    cut_flags.append(chg)
                Kratos.Logger.PrintInfo(
                    "CapDebug",
                    f"elem {elem.Id} face_nodes {face_nodes} node_ids {[n.Id for n in nodes]} "
                    f"phi_raw {[self.control_field[n.Id-1] for n in nodes]} "
                    f"phi_signed {F4s.tolist()} inside {inside_dbg.tolist()} cuts {cut_flags}"
                )
                self._caps_debug_count += 1

            triangles = self._cap_face_from_edges(P4, F4s, tol)
            for tri in triangles:
                # tri is (3,3) coords
                if self._tri_area(tri[0], tri[1], tri[2]) <= self.min_triangle_area:
                    continue
                node_ids = [get_or_create_node_id(pt) for pt in tri]
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)
                cond_id += 1
        return cond_id

