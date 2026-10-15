import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestJ2PlaneStrainThreeComponents(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KM.Model()
        self.mp = self.model.CreateModelPart('material')
        for i, (x, y) in enumerate(((0., 0.), (1., 0.), (0., 1.)), 1):
            self.mp.CreateNewNode(i, x, y, 0.)
        self.geometry = KM.Triangle2D3(*[self.mp.Nodes[i] for i in (1, 2, 3)])
        self.properties = self.mp.CreateNewProperties(1)
        for name, value in {'YOUNG_MODULUS': 67000., 'POISSON_RATIO': .0,
                            'YIELD_STRESS': 201., 'ISOTROPIC_HARDENING_MODULUS': 0.,
                            'EXPONENTIAL_SATURATION_YIELD_STRESS': 201.,
                            'HARDENING_EXPONENT': 0.}.items():
            self.properties[KM.KratosGlobals.GetVariable(name)] = value

    def make_law(self, reference=False):
        if reference:
            law = CLA.SmallStrainJ2PlasticityPlaneStrain2DLaw()
        else:
            law = CLA.SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw()
        law.InitializeMaterial(self.properties, self.geometry, KM.Vector([1./3]*3))
        return law

    def response(self, law, strain, commit=False):
        values = KM.ConstitutiveLawParameters()
        values.SetMaterialProperties(self.properties)
        values.SetElementGeometry(self.geometry)
        values.SetProcessInfo(self.mp.ProcessInfo)
        strain_vector = KM.Vector(strain)
        stress = KM.Vector(len(strain))
        tangent = KM.Matrix(len(strain), len(strain))
        values.SetStrainVector(strain_vector)
        values.SetStressVector(stress)
        values.SetConstitutiveMatrix(tangent)
        flags = KM.Flags()
        for flag in (KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN,
                     KM.ConstitutiveLaw.COMPUTE_STRESS,
                     KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR):
            flags.Set(flag, True)
        values.SetOptions(flags)
        law.CalculateMaterialResponseCauchy(values)
        result = (KM.Vector(stress), KM.Matrix(tangent))
        self.assertVectorAlmostEqual(strain_vector, KM.Vector(strain), 12)
        if commit:
            law.FinalizeMaterialResponseCauchy(values)
        return result

    def test_history_and_hardening_match_four_components(self):
        for hardening in (0., 1200.):
            self.properties[KM.KratosGlobals.GetVariable('ISOTROPIC_HARDENING_MODULUS')] = hardening
            law, reference = self.make_law(), self.make_law(True)
            for e in ([.001, -.0002, .0003], [.004, -.009, .003],
                      [.008, -.012, -.002], [.002, -.003, .001], [-.006, .004, -.003]):
                stress, tangent = self.response(law, e, True)
                stress4, tangent4 = self.response(reference, [e[0], e[1], 0., e[2]], True)
                for i, ii in enumerate((0, 1, 3)):
                    self.assertAlmostEqual(stress[i], stress4[ii], 10)
                    for j, jj in enumerate((0, 1, 3)):
                        self.assertAlmostEqual(tangent[i, j], tangent4[ii, jj], 8)
                self.assertAlmostEqual(law.GetValue(KM.KratosGlobals.GetVariable("ACCUMULATED_PLASTIC_STRAIN")),
                                       reference.GetValue(KM.KratosGlobals.GetVariable("ACCUMULATED_PLASTIC_STRAIN")), 12)
                self.assertAlmostEqual(law.GetValue(CLA.EQUIVALENT_PLASTIC_STRAIN),
                                       reference.GetValue(KM.KratosGlobals.GetVariable("ACCUMULATED_PLASTIC_STRAIN")), 12)

    def test_tangent_trial_state_and_clone(self):
        law = self.make_law()
        self.response(law, [.004, -.009, .003], True)
        old_state = law.GetValue(CLA.EQUIVALENT_PLASTIC_STRAIN)
        e = [.008, -.012, -.002]
        stress, tangent = self.response(law, e)
        for j in range(3):
            plus, minus = list(e), list(e)
            plus[j] += 1e-8
            minus[j] -= 1e-8
            sp, _ = self.response(law, plus)
            sm, _ = self.response(law, minus)
            for i in range(3):
                self.assertAlmostEqual(tangent[i, j], (sp[i]-sm[i])/2e-8, delta=2e-4)
        self.assertEqual(law.GetValue(CLA.EQUIVALENT_PLASTIC_STRAIN), old_state)
        clone = law.Clone()
        self.assertVectorAlmostEqual(self.response(clone, e)[0], stress, 10)
        self.response(clone, e, True)
        self.assertEqual(law.GetValue(CLA.EQUIVALENT_PLASTIC_STRAIN), old_state)
        serializer = KM.StreamSerializer()
        self.properties[KM.CONSTITUTIVE_LAW] = clone
        serializer.Save('material', self.mp)
        restored_model = KM.Model()
        restored_mp = restored_model.CreateModelPart('material')
        serializer.Load('material', restored_mp)
        restored = restored_mp.GetProperties()[1][KM.CONSTITUTIVE_LAW]
        self.assertVectorAlmostEqual(self.response(restored, [.002, -.003, .001])[0],
                                     self.response(clone, [.002, -.003, .001])[0], 10)

    def test_analytic_plane_strain_response(self):
        law = self.make_law()
        self.assertEqual(law.GetStrainSize(), 3)
        stress, _ = self.response(law, [.01, 0., 0.])
        self.assertVectorAlmostEqual(stress, KM.Vector([357.3333333333333, 156.3333333333333, 0.]), 10)


if __name__ == '__main__':
    KratosUnittest.main()
