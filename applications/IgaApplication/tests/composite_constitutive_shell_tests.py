#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application

# - ABD checks impose generalized-strain states and recover the reconstructed A/B/D (+ transverse shear) matrix and compare to the analytic reference
# - Pagano simply-supported plate central-deflection benchmark (asymmetric [0/90] and symmetric [0/90/0]) vs CLPT (for 3p shell) andFSDT-Whitney (for 6p shell) Navier

import json
import math
import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest

GV = KM.KratosGlobals.GetVariable
GLS = GV("GREEN_LAGRANGE_STRAIN_VECTOR")
MF = GV("MEMBRANE_FORCE")
IM = GV("INTERNAL_MOMENT")
SF1 = GV("SHEAR_FORCE_1")
SF2 = GV("SHEAR_FORCE_2")

# mechanics helpers 
def make_shell_model_part(model, name="IgaModelPart"):
    mp = model.CreateModelPart(name)
    mp.SetBufferSize(2)
    for var in (KM.DISPLACEMENT, KM.REACTION, KM.ROTATION, KM.REACTION_MOMENT,
                IGA.POINT_LOAD, IGA.DEAD_LOAD, SMA.POINT_LOAD):
        mp.AddNodalSolutionStepVariable(var)
    return mp

def _open_knot_vector_full(p, n_elem, length):
    step = length / float(n_elem)
    interior = [step * i for i in range(1, n_elem)]
    return [0.0] * (p + 1) + interior + [float(length)] * (p + 1)

def _greville(knots_full, p, n_cp, scale):
    span = knots_full[-1]
    return [scale * (sum(knots_full[i + 1:i + p + 1]) / float(p)) / span for i in range(n_cp)]

def make_flat_plate(model_part, A, B, p_u, p_v, n_elem_u, n_elem_v, shear=0.0):
    n_u, n_v = p_u + n_elem_u, p_v + n_elem_v
    Uf = _open_knot_vector_full(p_u, n_elem_u, A)
    Vf = _open_knot_vector_full(p_v, n_elem_v, B)
    xs = _greville(Uf, p_u, n_u, A)
    ys = _greville(Vf, p_v, n_v, B)
    nodes = KM.NodesVector()
    nid = max((n.Id for n in model_part.Nodes), default=0) + 1
    for vy in ys:
        for ux in xs:
            nodes.append(model_part.CreateNewNode(nid, ux + shear * vy, vy, 0.0))
            nid += 1
    ku = KM.Vector(len(Uf) - 2)
    for i, v in enumerate(Uf[1:-1]):
        ku[i] = v
    kv = KM.Vector(len(Vf) - 2)
    for i, v in enumerate(Vf[1:-1]):
        kv[i] = v
    return KM.NurbsSurfaceGeometry3D(nodes, p_u, p_v, ku, kv), n_u, n_v

def create_shell_elements(model_part, surface, element_name, props, deriv_order=3, eid_start=1):
    qpg = KM.GeometriesVector()
    surface.CreateQuadraturePointGeometries(qpg, deriv_order)
    for i in range(len(qpg)):
        model_part.CreateNewElement(element_name, eid_start + i, qpg[i], props)
    return qpg

def attach_composite_cl(props, angles_deg, thicknesses, plies, drill_E, drill_nu, density=1.0,
                        drilling_penalty=0.0):
    t_total = sum(thicknesses)
    z, z_lo = [], -0.5 * t_total
    for tk in thicknesses:
        z.append(z_lo + 0.5 * tk)
        z_lo += tk
    params = KM.Parameters(json.dumps({
        "z_layer_coordinate_vector": z,
        "Euler_angle_layer_vector": list(angles_deg),
        "thickness_layer_vector": list(thicknesses),
    }))
    props.SetValue(KM.THICKNESS, t_total)
    props.SetValue(KM.DENSITY, density)
    props.SetValue(KM.YOUNG_MODULUS, drill_E)
    props.SetValue(KM.POISSON_RATIO, drill_nu)
    if drilling_penalty != 0.0:
        props.SetValue(GV("DRILLING_PENALTY"), float(drilling_penalty))
    props.SetValue(KM.CONSTITUTIVE_LAW, IGA.IgaThicknessIntegratedCompositeLaw().Create(params))
    for k, (tk, (E1, E2, nu12, G12, G13, G23)) in enumerate(zip(thicknesses, plies)):
        sub = KM.Properties(props.Id * 100 + (k + 1))
        sub.SetValue(KM.DENSITY, density)
        sub.SetValue(KM.THICKNESS, tk)
        sub.SetValue(GV("YOUNG_MODULUS_X"), E1)
        sub.SetValue(GV("YOUNG_MODULUS_Y"), E2)
        sub.SetValue(GV("POISSON_RATIO_XY"), nu12)
        sub.SetValue(GV("SHEAR_MODULUS_XY"), G12)
        ortho = KM.Vector(6)
        ortho[0] = E1; ortho[1] = E2; ortho[2] = E2
        ortho[3] = nu12; ortho[4] = nu12; ortho[5] = nu12
        sub.SetValue(GV("ORTHOTROPIC_ELASTIC_CONSTANTS"), ortho)
        sub.SetValue(GV("SHEAR_MODULUS_XZ"), G13)
        sub.SetValue(GV("SHEAR_MODULUS_YZ"), G23)
        sub.SetValue(KM.CONSTITUTIVE_LAW, CLA.LinearElasticOrthotropic2DLaw())
        props.AddSubProperties(sub)


def apply_sinusoidal_load(main_mp, qp_geoms, P0, A, B, props_id=3, submp="SinusoidalLoad"):
    """q(x,y) = -P0 sin(pi x/A) sin(pi y/B) (downward) via one LoadCondition per QP geometry."""
    sub = main_mp.CreateSubModelPart(submp) if not main_mp.HasSubModelPart(submp) else main_mp.GetSubModelPart(submp)
    if not main_mp.HasProperties(props_id):
        main_mp.CreateNewProperties(props_id)
    props = main_mp.GetProperties()[props_id]
    cid = max((c.Id for c in main_mp.Conditions), default=0) + 1
    for i in range(len(qp_geoms)):
        c = qp_geoms[i].Center()
        q_load = -P0 * math.sin(math.pi * float(c[0]) / A) * math.sin(math.pi * float(c[1]) / B)
        cond = main_mp.CreateNewCondition("LoadCondition", cid + i, qp_geoms[i], props)
        f = KM.Vector(3); f[0] = 0.0; f[1] = 0.0; f[2] = q_load
        cond.SetValue(IGA.DEAD_LOAD, f)
        sub.AddCondition(cond)


def add_dofs(main_mp, add_rotation_dofs=True):
    pairs = [(KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
             (KM.DISPLACEMENT_Z, KM.REACTION_Z)]
    if add_rotation_dofs:
        pairs += [(KM.ROTATION_X, KM.REACTION_MOMENT_X), (KM.ROTATION_Y, KM.REACTION_MOMENT_Y),
                  (KM.ROTATION_Z, KM.REACTION_MOMENT_Z)]
    for d, r in pairs:
        KM.VariableUtils().AddDof(d, r, main_mp)


def solve_linear(main_mp, add_rotation_dofs=True):
    add_dofs(main_mp, add_rotation_dofs)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    linear_solver = KM.SkylineLUFactorizationSolver()
    bs = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
    bs.SetCalculateReactionsFlag(False)
    strat = KM.ResidualBasedLinearStrategy(main_mp, scheme, bs, False, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    main_mp.CloneTimeStep(1.0)
    strat.Solve()
    strat.FinalizeSolutionStep()
    return bool(strat.IsConverged())


def w_at_param(model, surface, p_u, p_v, knots_u, knots_v, nodes_of_surface):
    """IGA-interpolated DISPLACEMENT_Z at the patch parametric centre (build a w-surface
    whose CPs are (0,0,w_a) and evaluate its GlobalCoordinates at the centre)."""
    tmp = model.CreateModelPart("__wsurf_%d" % id(surface))
    wnodes = KM.NodesVector()
    nid = 1
    for nd in nodes_of_surface:
        wnodes.append(tmp.CreateNewNode(nid, 0.0, 0.0, float(nd.GetSolutionStepValue(KM.DISPLACEMENT_Z))))
        nid += 1
    ku = KM.Vector(len(knots_u))
    for i, v in enumerate(knots_u):
        ku[i] = float(v)
    kv = KM.Vector(len(knots_v))
    for i, v in enumerate(knots_v):
        kv[i] = float(v)
    wsurf = KM.NurbsSurfaceGeometry3D(wnodes, p_u, p_v, ku, kv)
    loc = KM.Array3()
    loc[0] = 0.5 * (knots_u[0] + knots_u[-1])
    loc[1] = 0.5 * (knots_v[0] + knots_v[-1])
    loc[2] = 0.0
    z = wsurf.GlobalCoordinates(loc)[2]
    model.DeleteModelPart("__wsurf_%d" % id(surface))
    return z


A_LEN = 1.0
H = A_LEN / 20.0
E2 = 7.0e9
PLY = (25.0 * E2, E2, 0.25, 0.5 * E2, 0.5 * E2, 0.2 * E2)
ANGLES = [0.0, 90.0]
FRACS = [0.5, 0.5]
P_DEG = 3
N_KNOTS = 4
C = 1.0e-3

REF_A = [[4561403508.77193, 87719298.24561404, 0.0],
         [87719298.24561404, 4561403508.77193, 0.0],
         [0.0, 0.0, 175000000.0]]
REF_B = [[-52631578.94736843, 0.0, 0.0],
         [0.0, 52631578.94736843, 0.0],
         [0.0, 0.0, 0.0]]
REF_D = [[950292.3976608189, 18274.853801169593, 0.0],
         [18274.853801169593, 950292.3976608189, 0.0],
         [0.0, 0.0, 36458.33333333334]]
REF_SH = [[100565191.67688099, 0.0], [0.0, 100565191.67688091]]

RTOL_ABD = 1e-8
RTOL_SHEAR = 1e-6


def _abd6_ref():
    A = np.array(REF_A); B = np.array(REF_B); D = np.array(REF_D)
    return np.block([[A, B], [B.T, D]])

FIELDS_6DOF = [
    lambda X, Y: ([C * X, 0, 0], [0, 0, 0]),
    lambda X, Y: ([0, C * Y, 0], [0, 0, 0]),
    lambda X, Y: ([0.5 * C * Y, 0.5 * C * X, 0], [0, 0, 0]),
    lambda X, Y: ([0, 0, 0], [0, C * X, 0]),
    lambda X, Y: ([0, 0, 0], [-C * Y, 0, 0]),
    lambda X, Y: ([0, 0, 0], [-C * X, C * Y, 0]),
    lambda X, Y: ([0, 0, C * X], [0, 0, 0]),
    lambda X, Y: ([0, 0, C * Y], [0, 0, 0]),
]

FIELDS_3P = [
    ("disp", lambda X, Y: (C * X, 0.0, 0.0)),
    ("disp", lambda X, Y: (0.0, C * Y, 0.0)),
    ("disp", lambda X, Y: (0.5 * C * Y, 0.5 * C * X, 0.0)),
    ("disp", lambda X, Y: (0.0, 0.0, 0.5 * C * X * X)),      
    ("disp", lambda X, Y: (0.0, 0.0, 0.5 * C * Y * Y)),       
    ("disp", lambda X, Y: (0.0, 0.0, C * X * Y)),             # chi_xy
]


def _build_patch(element_name):
    model = KM.Model()
    mp = make_shell_model_part(model)
    surf, n_u, n_v = make_flat_plate(mp, A_LEN, A_LEN, P_DEG, P_DEG, N_KNOTS, N_KNOTS)
    props = mp.GetProperties()[1]
    attach_composite_cl(props, ANGLES, [f * H for f in FRACS], [PLY, PLY],
                        drill_E=PLY[0], drill_nu=PLY[2], drilling_penalty=1.0e-5)
    is_kl = (element_name == "Shell3pElement")
    create_shell_elements(mp, surf, element_name, props, deriv_order=4 if is_kl else 3)
    info = mp.ProcessInfo
    for e in mp.Elements:
        e.Initialize(info)
    return mp, info, next(iter(mp.Elements)), [surf[i] for i in range(len(surf))]


def _recover_8(elem, info):
    eps = np.array(elem.CalculateOnIntegrationPoints(GLS, info)[0])
    N = np.array(elem.CalculateOnIntegrationPoints(MF, info)[0])
    M = np.array(elem.CalculateOnIntegrationPoints(IM, info)[0])
    Q = np.array([elem.CalculateOnIntegrationPoints(SF1, info)[0],
                  elem.CalculateOnIntegrationPoints(SF2, info)[0]])
    return np.hstack([eps[0:3], eps[3:6], eps[6:8]]), np.hstack([N, M, Q])


def _recover_6(elem, info):
    eps = np.array(elem.CalculateOnIntegrationPoints(GLS, info)[0])
    N = np.array(elem.CalculateOnIntegrationPoints(MF, info)[0])
    M = np.array(elem.CalculateOnIntegrationPoints(IM, info)[0])
    return np.hstack([eps[0:3], eps[3:6]]), np.hstack([N, M])


# Pagano simply-supported plate deflection benchmark
P0 = 1.0
PAGANO_P_DEG = 3
PAGANO_N_KNOTS = 12

PAGANO_LAMINATES = {
    "symmetric":  (1.0, 3.0, [0.0, 90.0, 0.0], [1 / 3., 1 / 3., 1 / 3.], [PLY, PLY, PLY]),
    "asymmetric": (1.0, 1.0, [0.0, 90.0],      [0.5, 0.5],               [PLY, PLY]),
}
CLPT_WBAR = {
    "symmetric":  {10: 0.503382, 20: 0.503382, 100: 0.503382},
    "asymmetric": {10: 1.063577, 20: 1.063577, 100: 1.063577},
}
WHITNEY_WBAR = {
    "symmetric":  {10: 0.931404, 20: 0.610577, 100: 0.507672},
    "asymmetric": {10: 1.239893, 20: 1.107656, 100: 1.06534},
}
RTOL_CLPT = 1.0e-2       # 3pvs CLPT Navier
RTOL_WHITNEY = 1.5e-2    # 6p vs FSDT/Whitney Navier
RTOL_VS_SOMMERWERK = 1.0e-4  


def _solve_wbar(element_name, lam_name, d_inv):
    """Solve one Pagano plate; return normalized centre deflection w_bar."""
    A, B, angles, fracs, plies = PAGANO_LAMINATES[lam_name]
    h = A / d_inv
    model = KM.Model()
    mp = make_shell_model_part(model)
    surf, n_u, n_v = make_flat_plate(mp, A, B, PAGANO_P_DEG, PAGANO_P_DEG, PAGANO_N_KNOTS, PAGANO_N_KNOTS)
    props = mp.GetProperties()[1]
    attach_composite_cl(props, angles, [f * h for f in fracs], plies,
                        drill_E=PLY[0], drill_nu=PLY[2], drilling_penalty=1.0e-5)
    is_kl = (element_name == "Shell3pElement")
    qpg = create_shell_elements(mp, surf, element_name, props, deriv_order=4 if is_kl else 3)

    # Hard simply-supported (SS-2) BCs matching the cross-ply Navier solution.
    def nd(v, uu):
        return surf[v * n_u + uu]
    for v in range(n_v):
        for uu in (0, n_u - 1):       # x = 0, A edges
            n = nd(v, uu)
            n.Fix(KM.DISPLACEMENT_Y); n.Fix(KM.DISPLACEMENT_Z)
            if not is_kl:
                n.Fix(KM.ROTATION_X); n.Fix(KM.ROTATION_Z)
    for uu in range(n_u):
        for v in (0, n_v - 1):        # y = 0, B edges
            n = nd(v, uu)
            n.Fix(KM.DISPLACEMENT_X); n.Fix(KM.DISPLACEMENT_Z)
            if not is_kl:
                n.Fix(KM.ROTATION_Y); n.Fix(KM.ROTATION_Z)

    apply_sinusoidal_load(mp, qpg, P0, A, B)
    solve_linear(mp, add_rotation_dofs=not is_kl)

    ku = [surf.KnotsU()[i] for i in range(surf.NumberOfKnotsU())]
    kv = [surf.KnotsV()[i] for i in range(surf.NumberOfKnotsV())]
    nodes_of_surface = [surf[i] for i in range(len(surf))]
    w_center = w_at_param(model, surf, PAGANO_P_DEG, PAGANO_P_DEG, ku, kv, nodes_of_surface)
    w_max = -w_center
    return 100.0 * h ** 3 * E2 / (P0 * A ** 4) * w_max


SKEW_E, SKEW_NU, SKEW_T = 1.0e7, 0.3, 0.01
SKEW_C = 1.0e-3


def _membrane_energy_ratio(element_name, shear):
    model = KM.Model()
    mp = make_shell_model_part(model)
    surf, n_u, n_v = make_flat_plate(mp, A_LEN, A_LEN, P_DEG, P_DEG, N_KNOTS, N_KNOTS, shear=shear)
    props = mp.GetProperties()[1]
    props.SetValue(KM.THICKNESS, SKEW_T)
    props.SetValue(KM.YOUNG_MODULUS, SKEW_E)
    props.SetValue(KM.POISSON_RATIO, SKEW_NU)
    props.SetValue(KM.DENSITY, 1.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStress2DLaw())
    create_shell_elements(mp, surf, element_name, props, deriv_order=3)
    add_dofs(mp, add_rotation_dofs=True)
    for i in range(len(surf)):
        n = surf[i]
        dx = SKEW_C * (n.X0 + 0.5 * n.Y0)
        dy = SKEW_C * (0.5 * n.X0 + n.Y0)
        n.SetSolutionStepValue(KM.DISPLACEMENT, [dx, dy, 0.0])
        n.Fix(KM.DISPLACEMENT_X); n.Fix(KM.DISPLACEMENT_Y); n.Fix(KM.DISPLACEMENT_Z)
        for rd in (KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z):
            n.SetSolutionStepValue(rd, 0.0); n.Fix(rd)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    bs = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
    strat = KM.ResidualBasedLinearStrategy(mp, scheme, bs, True, False, False, False)
    strat.SetEchoLevel(0); strat.Initialize(); mp.CloneTimeStep(1.0); strat.Solve()
    U_fe = 0.0
    for i in range(len(surf)):
        n = surf[i]
        d = n.GetSolutionStepValue(KM.DISPLACEMENT)
        r = n.GetSolutionStepValue(KM.REACTION)
        U_fe += 0.5 * (d[0] * r[0] + d[1] * r[1] + d[2] * r[2])
    fac = SKEW_E / (1.0 - SKEW_NU ** 2)
    eps_C_eps = fac * SKEW_C ** 2 * (2.0 + 2.0 * SKEW_NU + 0.5 * (1.0 - SKEW_NU))
    U_analytic = 0.5 * eps_C_eps * SKEW_T * (A_LEN * A_LEN)
    return U_fe / U_analytic


# tests
class CompositeConstitutiveShellTests(KratosUnittest.TestCase):

    def _check_6dof(self, element_name):
        mp, info, elem, nodes = _build_patch(element_name)
        Gm = np.zeros((8, 8)); Sm = np.zeros((8, 8))
        for k, field in enumerate(FIELDS_6DOF):
            for n in nodes:
                d, r = field(n.X0, n.Y0)
                n.SetSolutionStepValue(KM.DISPLACEMENT, [float(d[0]), float(d[1]), float(d[2])])
                n.SetSolutionStepValue(KM.ROTATION, [float(r[0]), float(r[1]), float(r[2])])
            g, s = _recover_8(elem, info)
            Gm[:, k] = g; Sm[:, k] = s
        C8 = Sm @ np.linalg.inv(Gm)
        abd_ref = _abd6_ref()
        err_abd = np.max(np.abs(C8[0:6, 0:6] - abd_ref)) / max(np.max(np.abs(abd_ref)), 1.0)
        sh_ref = np.array(REF_SH)
        err_sh = np.max(np.abs(C8[6:8, 6:8] - sh_ref)) / max(np.max(np.abs(sh_ref)), 1.0)
        self.assertLess(err_abd, RTOL_ABD, msg=f"{element_name}: ABD6 rel err {err_abd:.2e}")
        self.assertLess(err_sh, RTOL_SHEAR, msg=f"{element_name}: shear block rel err {err_sh:.2e}")

    def test_Shell6pElement_ABD(self):
        self._check_6dof("Shell6pElement")

    def test_Shell6pElementSommerwerk_ABD(self):
        self._check_6dof("Shell6pElement_Sommerwerk")

    def test_Shell6pSommerwerk_membrane_energy_orthogonal(self):
        """Baseline: constant-strain membrane energy is exact on an orthogonal
        parametrisation (a1 . a2 = 0)."""
        ratio = _membrane_energy_ratio("Shell6pElement_Sommerwerk", shear=0.0)
        self.assertAlmostEqual(ratio, 1.0, delta=1e-6,
            msg=f"orthogonal membrane-energy ratio {ratio:.6f} != 1")

    def test_Shell6pSommerwerk_membrane_energy_skewed(self):
        """On a skewed parametrisation (a1 . a2 != 0) the constant-strain membrane energy
        must still be exact. The element builds the membrane stiffness directly from the
        contravariant metric a^{ab} (Gl. 2.44) and uses an orthonormal local frame
        e_2 = e_3 x e_1 (= a^2/||a^2||). A diagonal length-normalised frame e_2 = a_2/||a_2||
        would ignore the skew and over-stiffen the membrane (ratio ~1.67 at shear=0.5)."""
        ratio = _membrane_energy_ratio("Shell6pElement_Sommerwerk", shear=0.5)
        self.assertAlmostEqual(ratio, 1.0, delta=1e-3,
            msg=f"skewed membrane-energy ratio {ratio:.6f} != 1")

    def test_Shell3pElement_ABD(self):
        """KL: no transverse shear -> recover the 6x6 [[A,B],[B^T,D]] only. The 3p reads node
        COORDINATES, so impose the states by moving control points."""
        mp, info, elem, nodes = _build_patch("Shell3pElement")
        Gm = np.zeros((6, 6)); Sm = np.zeros((6, 6))
        X0 = [n.X0 for n in nodes]; Y0 = [n.Y0 for n in nodes]
        for k, (_, field) in enumerate(FIELDS_3P):
            for i, n in enumerate(nodes):
                dx, dy, dz = field(X0[i], Y0[i])
                n.X = X0[i] + dx; n.Y = Y0[i] + dy; n.Z = dz
            g, s = _recover_6(elem, info)
            Gm[:, k] = g; Sm[:, k] = s
        C6 = Sm @ np.linalg.inv(Gm)
        abd_ref = _abd6_ref()
        err = np.max(np.abs(C6 - abd_ref)) / max(np.max(np.abs(abd_ref)), 1.0)
        self.assertLess(err, 1e-6, msg=f"Shell3pElement: ABD6 rel err {err:.2e}")

    def test_Shell3p_Pagano_CLPT(self):
        """KL Shell3pElement carries no transverse shear -> matches CLPT Navier."""
        for lam in ("symmetric", "asymmetric"):
            for d_inv, w_ref in CLPT_WBAR[lam].items():
                w_fe = _solve_wbar("Shell3pElement", lam, d_inv)
                self.assertAlmostEqual(w_fe / w_ref, 1.0, delta=RTOL_CLPT,
                    msg=f"3p {lam} d=1/{d_inv}: w_bar {w_fe:.4f} vs CLPT {w_ref:.4f}")

    def test_Shell6p_Pagano_Whitney_and_matches_Sommerwerk(self):
        """Legacy RM Shell6pElement matches FSDT/Whitney AND the Shell6pElement_Sommerwerk
        (both pull the same IgaThicknessIntegratedCompositeLaw)."""
        for lam in ("symmetric", "asymmetric"):
            for d_inv, w_ref in WHITNEY_WBAR[lam].items():
                w_leg = _solve_wbar("Shell6pElement", lam, d_inv)
                w_som = _solve_wbar("Shell6pElement_Sommerwerk", lam, d_inv)
                self.assertAlmostEqual(w_leg / w_ref, 1.0, delta=RTOL_WHITNEY,
                    msg=f"6p {lam} d=1/{d_inv}: w_bar {w_leg:.4f} vs Whitney {w_ref:.4f}")
                self.assertAlmostEqual(w_leg / w_som, 1.0, delta=RTOL_VS_SOMMERWERK,
                    msg=f"{lam} d=1/{d_inv}: legacy 6p {w_leg:.5f} != Sommerwerk {w_som:.5f}")


if __name__ == "__main__":
    KratosUnittest.main()
