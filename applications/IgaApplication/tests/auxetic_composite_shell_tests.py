import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from composite_constitutive_shell_tests import (
    GV, A_LEN, H, P_DEG, N_KNOTS, RTOL_ABD, FIELDS_6DOF,
    make_shell_model_part, make_flat_plate, attach_composite_cl,
    create_shell_elements, _recover_8,
)


def _build_single_layer(element_name, angles, thicknesses, plies, drill_nu=0.3):
    model = KM.Model()
    mp = make_shell_model_part(model)
    surf, n_u, n_v = make_flat_plate(mp, A_LEN, A_LEN, P_DEG, P_DEG, N_KNOTS, N_KNOTS)
    props = mp.GetProperties()[1]
    attach_composite_cl(props, angles, thicknesses, plies,
                        drill_E=plies[0][0], drill_nu=drill_nu, drilling_penalty=1.0e-5)
    is_kl = (element_name == "Shell3pElement")
    create_shell_elements(mp, surf, element_name, props, deriv_order=4 if is_kl else 3)
    info = mp.ProcessInfo
    for e in mp.Elements:
        e.Initialize(info)
    nodes = [surf[i] for i in range(len(surf))]
    return mp, info, next(iter(mp.Elements)), nodes


def _recover_C8(element_name, angles, thicknesses, plies):
    mp, info, elem, nodes = _build_single_layer(element_name, angles, thicknesses, plies)
    Gm = np.zeros((8, 8)); Sm = np.zeros((8, 8))
    for k, field in enumerate(FIELDS_6DOF):
        for n in nodes:
            d, r = field(n.X0, n.Y0)
            n.SetSolutionStepValue(KM.DISPLACEMENT, [float(d[0]), float(d[1]), float(d[2])])
            n.SetSolutionStepValue(KM.ROTATION, [float(r[0]), float(r[1]), float(r[2])])
        g, s = _recover_8(elem, info)
        Gm[:, k] = g; Sm[:, k] = s
    return Sm @ np.linalg.inv(Gm)


class AuxeticCompositeShellTests(KratosUnittest.TestCase):

    def test_Shell6p_negative_poisson_inplane_coupling(self):
        E0, nu, t = 70.0e9, -0.5, H
        G0 = E0 / (2.0 * (1.0 + nu))
        ply = (E0, E0, nu, G0, G0, G0)
        fac = E0 * t / (1.0 - nu ** 2)
        A_ref = np.array([[fac, fac * nu, 0.0],
                          [fac * nu, fac, 0.0],
                          [0.0, 0.0, 0.5 * (1.0 - nu) * fac]])
        for element_name in ("Shell6pElement", "Shell6pElement_Sommerwerk"):
            C8 = _recover_C8(element_name, [0.0], [t], [ply])
            A = C8[0:3, 0:3]
            self.assertLess(A[0, 1], 0.0,
                msg=f"{element_name}: expected negative Poisson coupling A12, got {A[0, 1]:.3e}")
            err = np.max(np.abs(A - A_ref)) / np.max(np.abs(A_ref))
            self.assertLess(err, RTOL_ABD,
                msg=f"{element_name}: auxetic membrane A rel err {err:.2e}")

    def test_orthotropic_ply_stability_check(self):
        mp, info, elem, _ = _build_single_layer(
            "Shell6pElement", [0.0], [H], [(70.0e9, 70.0e9, -0.5, 23.0e9, 1.0, 1.0)])
        geom = elem.GetGeometry()
        law = CLA.LinearElasticOrthotropic2DLaw()

        def props(E1, E2, nu12, G12):
            p = KM.Properties(99)
            p.SetValue(GV("YOUNG_MODULUS_X"), E1)
            p.SetValue(GV("YOUNG_MODULUS_Y"), E2)
            p.SetValue(GV("POISSON_RATIO_XY"), nu12)
            p.SetValue(GV("SHEAR_MODULUS_XY"), G12)
            p.SetValue(KM.DENSITY, 1.0)
            return p

        self.assertEqual(law.Check(props(70e9, 70e9, -0.5, 23e9), geom, info), 0)
        with self.assertRaises(RuntimeError):
            law.Check(props(1.0e9, 100.0e9, 0.5, 1.0e9), geom, info)
        with self.assertRaises(RuntimeError):
            law.Check(props(70e9, 70e9, -1.2, 1.0e9), geom, info)
        with self.assertRaises(RuntimeError):
            law.Check(props(70e9, 70e9, -0.3, -1.0), geom, info)


if __name__ == "__main__":
    KratosUnittest.main()
