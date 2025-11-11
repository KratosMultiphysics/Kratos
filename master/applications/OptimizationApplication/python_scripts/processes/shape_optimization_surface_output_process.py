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
                "zero_tol"                    : 1e-12
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
            key = tuple(np.round(coord, 9))
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
                node_ids = [get_or_create_node_id(verts_phys[i]) for i in tri]
                mp_interface.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)
                cond_id += 1
            
        # Add caps on domain boundary faces so the surface closes
        if self.cap_domain_boundary:
            cond_id = self._add_boundary_caps_hex8(mp_interface, get_or_create_node_id, cond_id, prop)

        Kratos.Logger.PrintInfo("InterfaceOutputProcess", f"✅ Created {len(self.created_coordinates)} interface nodes and {cond_id-1} surface conditions.")

    def _interp_on_edge(self, p0, p1, f0, f1):
        # robust linear interpolation for φ=0 along an edge
        denom = (f0 - f1)
        if abs(denom) < 1e-30:
            t = 0.5
        else:
            t = f0 / denom
        return p0 + t * (p1 - p0)

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

    def _add_boundary_caps_hex8(self, mp_interface, get_or_create_node_id, cond_id, prop):
        """
        For each domain boundary face, add triangles covering the subset where φ<=0 (or φ>=0)
        to produce watertight closure.
        """
        sgn = 1.0 if self.solid_is_negative else -1.0   # treat inside as φ*sgn <= 0
        TRI_IN_FACE = [(0,1,2), (0,2,3)]  # split quad to 2 tris

        for elem, face_nodes in self._collect_boundary_faces_hex8():
            geom = elem.GetGeometry()

            # gather face data in local face order
            nodes = [geom[i] for i in face_nodes]
            P4 = np.array([[n.X, n.Y, n.Z] for n in nodes])
            F4 = np.array([sgn * self.control_field[n.Id - 1] for n in nodes])

            # clip each of the two triangles by φ<=0 and emit triangles
            for (a,b,c) in TRI_IN_FACE:
                Ptri = P4[[a,b,c], :]
                Ftri = F4[[a,b,c]]
                clipped = self._clip_triangle_phi_le_zero(Ptri, Ftri)
                for tri_pts in clipped:
                    node_ids = [get_or_create_node_id(pt) for pt in tri_pts]
                    mp_interface.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)
                    cond_id += 1

        return cond_id

