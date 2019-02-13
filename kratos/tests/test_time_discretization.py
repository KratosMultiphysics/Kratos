from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTimeDiscretization(KratosUnittest.TestCase):

    def _test_bdf_time_discretization(self, bdf, exp_results, *delta_times):
        test_model = KM.Model()
        model_part = test_model.CreateModelPart("test_mp")
        for dt in reversed(delta_times):
            model_part.CloneTimeStep(model_part.ProcessInfo[KM.TIME] + dt) # filling the time-step-info

        exp_size = len(exp_results)

        coeffs = bdf.ComputeBDFCoefficients(*delta_times)

        self.assertEqual(len(coeffs), exp_size)
        for coeff, exp_coeff in zip(coeffs, exp_results):
            self.assertAlmostEqual(coeff, exp_coeff)

        self.assertEqual(KM.GetMinimumBufferSize(bdf), exp_size)

        # test with passing a ProcessInfo
        coeffs = bdf.ComputeBDFCoefficients(model_part.ProcessInfo)
        self.assertEqual(len(coeffs), exp_size)
        for coeff, exp_coeff in zip(coeffs, exp_results):
            self.assertAlmostEqual(coeff, exp_coeff)

    def test_BDF1(self):
        bdf = KM.BDF1()

        delta_time = 0.11
        exp_results = [
            1.0/delta_time,
            -1.0/delta_time
        ]

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF2(self):
        bdf = KM.BDF2()

        delta_time = 0.11
        prev_delta_time = 0.089
        rho = prev_delta_time / delta_time;
        time_coeff = 1.0 / (delta_time * rho * rho + delta_time * rho);
        exp_results = [
            time_coeff * (rho * rho + 2.0 * rho),
            -time_coeff * (rho * rho + 2.0 * rho + 1.0),
            time_coeff
        ]

        self._test_bdf_time_discretization(bdf, exp_results, delta_time, prev_delta_time)

    def test_BDF3(self):
        bdf = KM.BDF3()

        delta_time = 0.11
        exp_results = [
            11.0/(6.0*delta_time),
            -18.0/(6.0*delta_time),
            9.0/(6.0*delta_time),
            -2.0/(6.0*delta_time)
        ]

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF4(self):
        bdf = KM.BDF4()

        delta_time = 0.11
        exp_results = [
            25.0/(12.0*delta_time),
            -48.0/(12.0*delta_time),
            36.0/(12.0*delta_time),
            -16.0/(12.0*delta_time),
            3.0/(12.0*delta_time)
        ]

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF5(self):
        bdf = KM.BDF5()

        delta_time = 0.11
        exp_results = [
            137.0/(60.0*delta_time),
            -300.0/(60.0*delta_time),
            300.0/(60.0*delta_time),
            -200.0/(60.0*delta_time),
            75.0/(60.0*delta_time),
            -12.0/(60.0*delta_time)
        ]

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF6(self):
        bdf = KM.BDF6()

        delta_time = 0.11
        exp_results = [
            147.0/(60.0*delta_time),
            -360.0/(60.0*delta_time),
            450.0/(60.0*delta_time),
            -400.0/(60.0*delta_time),
            225.0/(60.0*delta_time),
            -72.0/(60.0*delta_time),
            10.0/(60.0*delta_time)
        ]

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_Newmark(self):
        gen_alpha = KM.Newmark()
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.25)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.5)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

    def test_Bossak(self):
        gen_alpha = KM.Bossak()
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.4225)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.8)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

    def test_GeneralizedAlpha(self):
        gen_alpha = KM.GeneralizedAlpha()

        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.4225)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.8)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)
        self.assertAlmostEqual(gen_alpha.GetAlphaF(), 0.0)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

if __name__ == '__main__':
    KratosUnittest.main()
