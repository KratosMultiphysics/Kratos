import math

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    from KratosMultiphysics import ConstitutiveLawsApplication

@KratosUnittest.skipUnless(CheckIfApplicationsAvailable("ConstitutiveLawsApplication"), "ConstitutiveLawsApplication is not available")
class TestAlmansiStrainOutput(KratosUnittest.TestCase):
    """The Almansi strain output of the solid elements must be the Euler-Almansi strain e = 0.5 * (I - inv(F * F^T)), not the Green-Lagrange one"""

    # The system is not used by the scheme here, so empty containers are enough
    _A = KratosMultiphysics.CompressedMatrix()
    _Dx = KratosMultiphysics.Vector()
    _b = KratosMultiphysics.Vector()

    ELEMENTS = ["TotalLagrangianElement3D8N", "UpdatedLagrangianElement3D8N"]

    def _create_model_part(self, element_name):
        self.model = KratosMultiphysics.Model()
        mp = self.model.CreateModelPart("Main", 2)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        coordinates = [
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
        for i, (x, y, z) in enumerate(coordinates):
            mp.CreateNewNode(i + 1, x, y, z)

        props = mp.GetProperties()[1]
        props.SetValue(KratosMultiphysics.YOUNG_MODULUS, 210000.0)
        props.SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        props.SetValue(KratosMultiphysics.DENSITY, 1.0)
        props.SetValue(KratosMultiphysics.VOLUME_ACCELERATION, [0.0, 0.0, 0.0])
        props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, ConstitutiveLawsApplication.HyperElastic3DLaw())

        mp.CreateNewElement(element_name, 1, [1, 2, 3, 4, 5, 6, 7, 8], props)

        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        mp.ProcessInfo[KratosMultiphysics.STEP] = 0
        mp.ProcessInfo[KratosMultiphysics.TIME] = 0.0

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        scheme.InitializeElements(mp)
        return mp, scheme

    def _apply_deformation(self, mp, scheme, deformation_gradient):
        """Applies the homogeneous displacement field u = (F - I) * X in the first step"""
        mp.CloneTimeStep(1.0)
        mp.ProcessInfo[KratosMultiphysics.STEP] = 1
        for node in mp.Nodes:
            X0 = [node.X0, node.Y0, node.Z0]
            displacement = [sum((deformation_gradient[i][j] - (1.0 if i == j else 0.0)) * X0[j] for j in range(3)) for i in range(3)]
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, displacement)
            # The elements work with the current coordinates
            node.X = X0[0] + displacement[0]
            node.Y = X0[1] + displacement[1]
            node.Z = X0[2] + displacement[2]
        scheme.InitializeSolutionStep(mp, self._A, self._Dx, self._b)

    @staticmethod
    def _rotated_axial_stretch(stretch, angle):
        """F = R * diag(stretch, 1, 1) and the Almansi strain tensor R * diag(e, 0, 0) * R^T"""
        c = math.cos(angle)
        s = math.sin(angle)
        deformation_gradient = [[c * stretch, -s, 0.0], [s * stretch, c, 0.0], [0.0, 0.0, 1.0]]
        e = 0.5 * (1.0 - 1.0 / stretch**2)
        almansi = [[c * c * e, c * s * e, 0.0], [c * s * e, s * s * e, 0.0], [0.0, 0.0, 0.0]]
        return deformation_gradient, almansi

    def _deformation_cases(self):
        """List of (name, deformation gradient, Almansi strain tensor)"""
        cases = []

        # Axial stretches
        for stretch, reference in [(1.002, 0.00199401596010), (1.004, 0.00397612736306)]:
            deformation_gradient = [[stretch, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            almansi = [[reference, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
            cases.append(("axial stretch {}".format(stretch), deformation_gradient, almansi))

        # Simple shear
        gamma = 0.15
        deformation_gradient = [[1.0, gamma, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        almansi = [[0.0, 0.5 * gamma, 0.0], [0.5 * gamma, -0.5 * gamma**2, 0.0], [0.0, 0.0, 0.0]]
        cases.append(("simple shear", deformation_gradient, almansi))

        # Axial stretch with a superposed rigid rotation
        deformation_gradient, almansi = self._rotated_axial_stretch(1.004, math.radians(30.0))
        cases.append(("rotated axial stretch", deformation_gradient, almansi))

        return cases

    @staticmethod
    def _to_voigt_strain(tensor):
        # xx, yy, zz, xy, yz, xz (engineering shear strains)
        return [tensor[0][0], tensor[1][1], tensor[2][2], 2.0 * tensor[0][1], 2.0 * tensor[1][2], 2.0 * tensor[0][2]]

    def test_AlmansiStrain(self):
        for element_name in self.ELEMENTS:
            for case_name, deformation_gradient, almansi in self._deformation_cases():
                with self.subTest(element=element_name, case=case_name):
                    mp, scheme = self._create_model_part(element_name)
                    self._apply_deformation(mp, scheme, deformation_gradient)
                    element = mp.GetElement(1)

                    reference_vector = self._to_voigt_strain(almansi)
                    vectors = element.CalculateOnIntegrationPoints(KratosMultiphysics.ALMANSI_STRAIN_VECTOR, mp.ProcessInfo)
                    tensors = element.CalculateOnIntegrationPoints(KratosMultiphysics.ALMANSI_STRAIN_TENSOR, mp.ProcessInfo)
                    self.assertEqual(len(vectors), 8)
                    self.assertEqual(len(tensors), 8)
                    for vector, tensor in zip(vectors, tensors):
                        self.assertVectorAlmostEqual(vector, reference_vector, 9)
                        for i in range(3):
                            for j in range(3):
                                self.assertAlmostEqual(tensor[i, j], almansi[i][j], 9)

    def test_GreenLagrangeStrain(self):
        """The Green-Lagrange output must remain E = 0.5 * (F^T * F - I)"""
        for element_name in self.ELEMENTS:
            for case_name, F, _ in self._deformation_cases():
                with self.subTest(element=element_name, case=case_name):
                    mp, scheme = self._create_model_part(element_name)
                    self._apply_deformation(mp, scheme, F)

                    green_lagrange = [[0.5 * (sum(F[k][i] * F[k][j] for k in range(3)) - (1.0 if i == j else 0.0)) for j in range(3)] for i in range(3)]
                    reference_vector = self._to_voigt_strain(green_lagrange)
                    vectors = mp.GetElement(1).CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, mp.ProcessInfo)
                    for vector in vectors:
                        self.assertVectorAlmostEqual(vector, reference_vector, 9)

    def test_RepeatedOutputRequests(self):
        """Repeated output requests must not alter the results"""
        for element_name in self.ELEMENTS:
            with self.subTest(element=element_name):
                _, F, almansi = self._deformation_cases()[-1]
                mp, scheme = self._create_model_part(element_name)
                self._apply_deformation(mp, scheme, F)
                reference_vector = self._to_voigt_strain(almansi)
                for _ in range(3):
                    for vector in mp.GetElement(1).CalculateOnIntegrationPoints(KratosMultiphysics.ALMANSI_STRAIN_VECTOR, mp.ProcessInfo):
                        self.assertVectorAlmostEqual(vector, reference_vector, 9)

if __name__ == '__main__':
    KratosUnittest.main()
