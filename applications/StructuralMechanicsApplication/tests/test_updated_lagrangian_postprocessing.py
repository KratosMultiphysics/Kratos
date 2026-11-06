import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    from KratosMultiphysics import ConstitutiveLawsApplication

@KratosUnittest.skipUnless(CheckIfApplicationsAvailable("ConstitutiveLawsApplication"), "ConstitutiveLawsApplication is not available")
class TestUpdatedLagrangianPostprocessing(KratosUnittest.TestCase):
    """Integration point output must not change after FinalizeSolutionStep if the displacements do not change"""

    # Imposed axial stretch at each step and the corresponding reference axial Cauchy stress
    # The system is not used by the scheme here, so empty containers are enough
    _A = KratosMultiphysics.CompressedMatrix()
    _Dx = KratosMultiphysics.Vector()
    _b = KratosMultiphysics.Vector()

    STRETCHES = [1.002, 1.004]
    REFERENCE_STRESS = [564.337033, 1126.588530]

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

    def _start_step(self, mp, scheme, step, stretch):
        mp.CloneTimeStep(float(step))
        mp.ProcessInfo[KratosMultiphysics.STEP] = step
        for node in mp.Nodes:
            displacement = (stretch - 1.0) * node.X0
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, displacement)
            # The elements work with the current coordinates
            node.X = node.X0 + displacement
        scheme.InitializeSolutionStep(mp, self._A, self._Dx, self._b)

    @staticmethod
    def _stresses(mp):
        return mp.GetElement(1).CalculateOnIntegrationPoints(KratosMultiphysics.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)

    def _check_against_reference(self, stresses, step_index):
        for stress in stresses:
            self.assertAlmostEqual(stress[0], self.REFERENCE_STRESS[step_index], 4)

    def _run(self, element_name):
        mp, scheme = self._create_model_part(element_name)
        for i, stretch in enumerate(self.STRETCHES):
            self._start_step(mp, scheme, i + 1, stretch)

            before = self._stresses(mp)
            self._check_against_reference(before, i)

            scheme.FinalizeSolutionStep(mp, self._A, self._Dx, self._b)

            # Several queries after finalization must be identical and equal to the ones before it
            for _ in range(2):
                after = self._stresses(mp)
                self._check_against_reference(after, i)
                for stress_before, stress_after in zip(before, after):
                    self.assertVectorAlmostEqual(stress_before, stress_after, 8)

    def test_TotalLagrangian(self):
        self._run("TotalLagrangianElement3D8N")

    def test_UpdatedLagrangian(self):
        self._run("UpdatedLagrangianElement3D8N")

if __name__ == '__main__':
    KratosUnittest.main()
