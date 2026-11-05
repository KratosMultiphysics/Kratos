"""Numerical checks of element-center von Mises stress recovery."""

from math import cos, sin, sqrt

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.StructuralMechanicsApplication.helper.responses import (
    get_von_mises,
)


class TestHelperResponses(KratosUnittest.TestCase):
    young = 210000.0
    poisson = 0.3
    triangle = ((0, 0, 0), (1, 0, 0), (0, 1, 0))
    quadrilateral = ((-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0))
    tetrahedron = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    hexahedron = ((-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                  (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1))

    def _make_part(self, name, coordinates, law, displacement, integration_order=None):
        model = KM.Model()
        part = model.CreateModelPart("Structure")
        part.SetBufferSize(2)
        part.ProcessInfo[KM.DOMAIN_SIZE] = 2 if "2D" in name else 3
        for variable in (KM.DISPLACEMENT, KM.ROTATION, KM.VOLUME_ACCELERATION):
            part.AddNodalSolutionStepVariable(variable)
        properties = part.CreateNewProperties(1)
        properties[KM.YOUNG_MODULUS] = self.young
        properties[KM.POISSON_RATIO] = self.poisson
        properties[KM.THICKNESS] = 0.1
        properties[KM.DENSITY] = 1.0
        properties[KM.CONSTITUTIVE_LAW] = law
        if integration_order is not None:
            properties[SMA.INTEGRATION_ORDER] = integration_order
        for index, coordinates_i in enumerate(coordinates, start=1):
            node = part.CreateNewNode(index, *coordinates_i)
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement(*coordinates_i))
            node.SetSolutionStepValue(KM.ROTATION, [0.0, 0.0, 0.0])
        element = part.CreateNewElement(
            name, 7, list(range(1, len(coordinates) + 1)), properties
        )
        element.Initialize(part.ProcessInfo)
        return model, part, element

    def _expected_plane_stress(self, xx, yy, engineering_xy):
        factor = self.young / (1.0 - self.poisson**2)
        stress_xx = factor * (xx + self.poisson * yy)
        stress_yy = factor * (yy + self.poisson * xx)
        stress_xy = self.young / (2.0 * (1.0 + self.poisson)) * engineering_xy
        return sqrt(stress_xx**2 + stress_yy**2 - stress_xx * stress_yy
                    + 3.0 * stress_xy**2)

    def test_plane_stress_triangle(self):
        _, part, _ = self._make_part(
            "SmallDisplacementElement2D3N", self.triangle,
            SMA.LinearElasticPlaneStress2DLaw(),
            lambda x, y, z: [0.001 * x + 0.0006 * y, -0.0002 * y, 0.0],
        )
        expected = self._expected_plane_stress(0.001, -0.0002, 0.0006)
        for modes in ((True, False, False), (False, True, False), (False, False, True)):
            with self.subTest(modes=modes):
                self.assertAlmostEqual(get_von_mises(*modes, part)[7], expected, 7)

    def test_plane_strain_includes_out_of_plane_stress(self):
        _, part, _ = self._make_part(
            "SmallDisplacementElement2D3N", self.triangle,
            SMA.LinearElasticPlaneStrain2DLaw(),
            lambda x, y, z: [0.001 * x, 0.001 * y, 0.0],
        )
        # Equibiaxial strain: sigma_xx = sigma_yy, sigma_zz = 2 nu sigma_xx.
        expected = self.young * 0.001 / (1.0 + self.poisson)
        self.assertAlmostEqual(get_von_mises(True, False, False, part)[7], expected, 7)

    def test_solid_tetrahedron_and_all_element_ids(self):
        _, part, element = self._make_part(
            "SmallDisplacementElement3D4N", self.tetrahedron,
            SMA.LinearElastic3DLaw(), lambda x, y, z: [0.002 * y, 0.0, 0.0],
        )
        other = part.CreateNewElement(
            "SmallDisplacementElement3D4N", 23, [1, 2, 3, 4], element.Properties
        )
        other.Initialize(part.ProcessInfo)
        values = get_von_mises(True, False, False, part)
        self.assertEqual(set(values), {7, 23})
        expected = sqrt(3.0) * self.young / (2.0 * (1.0 + self.poisson)) * 0.002
        for value in values.values():
            self.assertAlmostEqual(value, expected, 7)

    def test_hydrostatic_solid_has_zero_von_mises(self):
        _, part, _ = self._make_part(
            "SmallDisplacementElement3D8N", self.hexahedron,
            SMA.LinearElastic3DLaw(), lambda x, y, z: [0.001*x, 0.001*y, 0.001*z],
        )
        for modes in ((False, True, False), (False, False, True)):
            with self.subTest(modes=modes):
                self.assertAlmostEqual(get_von_mises(*modes, part)[7], 0.0, 8)

    def test_center_gp_uses_element_integration_order(self):
        _, part, _ = self._make_part(
            "SmallDisplacementElement3D8N", self.hexahedron,
            SMA.LinearElastic3DLaw(), lambda x, y, z: [0.001*y, 0.0, 0.0],
            integration_order=3,
        )
        expected = sqrt(3.0) * self.young / (2.0 * (1.0 + self.poisson)) * 0.001
        self.assertAlmostEqual(get_von_mises(True, False, False, part)[7], expected, 7)

    def test_prism_and_pyramid_solid_interpolation(self):
        prism = ((0, 0, -1), (1, 0, -1), (0, 1, -1),
                 (0, 0, 1), (1, 0, 1), (0, 1, 1))
        pyramid = ((-1, -1, -1), (1, -1, -1), (1, 1, -1),
                   (-1, 1, -1), (0, 0, 1))
        expected = sqrt(3.0) * self.young / (2.0 * (1.0 + self.poisson)) * 0.001
        for name, coordinates in (("SmallDisplacementElement3D6N", prism),
                                  ("SmallDisplacementElement3D5N", pyramid)):
            with self.subTest(element=name):
                _, part, _ = self._make_part(
                    name, coordinates, SMA.LinearElastic3DLaw(),
                    lambda x, y, z: [0.001 * y, 0.0, 0.0],
                )
                self.assertAlmostEqual(
                    get_von_mises(False, False, True, part)[7], expected, 7
                )

    def test_average_stress_before_von_mises(self):
        _, part, element = self._make_part(
            "SmallDisplacementElement2D4N", self.quadrilateral,
            SMA.LinearElasticPlaneStress2DLaw(),
            lambda x, y, z: [0.001 * x * y, 0.0, 0.0],
        )
        gp_vm = element.CalculateOnIntegrationPoints(
            SMA.VON_MISES_STRESS, part.ProcessInfo
        )
        self.assertGreater(sum(gp_vm) / len(gp_vm), 1.0)
        self.assertAlmostEqual(get_von_mises(False, True, False, part)[7], 0.0, 8)

    def test_interpolation_differs_from_gp_average(self):
        coordinates = self.quadrilateral + ((0, -1, 0), (1, 0, 0),
                                             (0, 1, 0), (-1, 0, 0))
        _, part, _ = self._make_part(
            "SmallDisplacementElement2D8N", coordinates,
            SMA.LinearElasticPlaneStress2DLaw(),
            lambda x, y, z: [0.001 * x * y**2, 0.0, 0.0], integration_order=4,
        )
        # The exact center strain is zero. The noncentral GPs sample nonzero y^2.
        self.assertAlmostEqual(get_von_mises(False, False, True, part)[7], 0.0, 6)
        self.assertGreater(get_von_mises(False, True, False, part)[7], 1.0)

    def test_interpolation_on_distorted_quadrilateral(self):
        coordinates = ((0, 0, 0), (3, 0.1, 0), (2.4, 2, 0), (-0.3, 1.2, 0))
        _, part, _ = self._make_part(
            "SmallDisplacementElement2D4N", coordinates,
            SMA.LinearElasticPlaneStress2DLaw(),
            lambda x, y, z: [0.001 * x, -0.0002 * y, 0.0],
        )
        expected = self._expected_plane_stress(0.001, -0.0002, 0.0)
        self.assertAlmostEqual(get_von_mises(False, False, True, part)[7], expected, 7)

    def test_missing_center_gp_is_explicit(self):
        _, part, _ = self._make_part(
            "SmallDisplacementElement2D4N", self.quadrilateral,
            SMA.LinearElasticPlaneStress2DLaw(), lambda x, y, z: [0.0, 0.0, 0.0],
        )
        with self.assertRaisesRegex(ValueError, "Element 7 .*no GP at its center"):
            get_von_mises(True, False, False, part)

    def test_linear_shell_middle_surface_in_rotated_plane(self):
        angle = 0.7
        expected = self._expected_plane_stress(0.001, 0.0, 0.0)
        for name, base in (("ShellThinElement3D3N", self.triangle),
                           ("LegacyThickShellElement3D4N", self.quadrilateral)):
            with self.subTest(element=name):
                coordinates = [(x, y*cos(angle), y*sin(angle)) for x, y, z in base]
                _, part, _ = self._make_part(
                    name, coordinates, SMA.LinearElasticPlaneStress2DLaw(),
                    lambda x, y, z: [0.001 * x, 0.0, 0.0],
                )
                for modes in ((False, True, False), (False, False, True)):
                    self.assertAlmostEqual(get_von_mises(*modes, part)[7], expected, 6)

    def test_thick_shell_middle_surface_includes_transverse_shear(self):
        _, part, element = self._make_part(
            "LegacyThickShellElement3D4N", self.quadrilateral,
            SMA.LinearElasticPlaneStress2DLaw(), lambda x, y, z: [0.0, 0.0, 0.001*x],
        )
        tensors = element.CalculateOnIntegrationPoints(
            SMA.SHELL_STRESS_MIDDLE_SURFACE_GLOBAL, part.ProcessInfo
        )
        shear_xz = sum(tensor[0, 2] for tensor in tensors) / len(tensors)
        shear_yz = sum(tensor[1, 2] for tensor in tensors) / len(tensors)
        shear_magnitude = sqrt(shear_xz**2 + shear_yz**2)
        self.assertGreater(shear_magnitude, 0.01)
        self.assertAlmostEqual(
            get_von_mises(False, True, False, part)[7], sqrt(3.0) * shear_magnitude, 7
        )

    def test_unsupported_elements_fail_before_stress_requests(self):
        for name, coordinates, message in (
                ("TrussElement3D2N", self.triangle[:2], "unsupported"),
                ("TotalLagrangianElement3D4N", self.tetrahedron, "unsupported"),
                ("IsotropicShellElement3D3N", self.triangle, "unsupported"),
                ("MITCThickShellElement3D4N", self.quadrilateral, "does not expose"),
                ("CSDSG3ThickShellLinearElement3D3N", self.triangle,
                 "does not expose")):
            with self.subTest(element=name):
                model = KM.Model()
                part = model.CreateModelPart("Unsupported")
                properties = part.CreateNewProperties(1)
                for index, coordinates_i in enumerate(coordinates, start=1):
                    part.CreateNewNode(index, *coordinates_i)
                part.CreateNewElement(
                    name, 7, list(range(1, len(coordinates) + 1)), properties
                )
                with self.assertRaisesRegex(NotImplementedError, message):
                    get_von_mises(False, True, False, part)

    def test_mode_selection_and_empty_model_part(self):
        model = KM.Model()
        part = model.CreateModelPart("Empty")
        self.assertEqual(get_von_mises(False, True, False, part), {})
        for modes in ((False, False, False), (True, True, False), (True, True, True)):
            with self.subTest(modes=modes):
                with self.assertRaisesRegex(ValueError, "exactly one"):
                    get_von_mises(*modes, part)
        with self.assertRaisesRegex(TypeError, "Booleans"):
            get_von_mises("true", False, False, part)


if __name__ == "__main__":
    KratosUnittest.main()
