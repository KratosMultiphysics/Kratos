import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTimeDiscretization(KratosUnittest.TestCase):
    @classmethod
    def _get_expected_results(self, time_order, delta_time, prev_delta_time = None):
        if time_order == 1:
            return [1.0/delta_time, -1.0/delta_time]
        elif time_order == 2:
            rho = prev_delta_time / delta_time;
            time_coeff = 1.0 / (delta_time * rho * rho + delta_time * rho);
            return [time_coeff * (rho * rho + 2.0 * rho), -time_coeff * (rho * rho + 2.0 * rho + 1.0), time_coeff]
        elif time_order == 3:
            return [11.0/(6.0*delta_time), -18.0/(6.0*delta_time), 9.0/(6.0*delta_time), -2.0/(6.0*delta_time)]
        elif time_order == 4:
            return [25.0/(12.0*delta_time), -48.0/(12.0*delta_time), 36.0/(12.0*delta_time), -16.0/(12.0*delta_time), 3.0/(12.0*delta_time)]
        elif time_order == 5:
            return [137.0/(60.0*delta_time), -300.0/(60.0*delta_time), 300.0/(60.0*delta_time), -200.0/(60.0*delta_time), 75.0/(60.0*delta_time), -12.0/(60.0*delta_time)]
        elif time_order == 6:
            return [147.0/(60.0*delta_time), -360.0/(60.0*delta_time), 450.0/(60.0*delta_time), -400.0/(60.0*delta_time), 225.0/(60.0*delta_time), -72.0/(60.0*delta_time), 10.0/(60.0*delta_time)]

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

        self.assertEqual(KM.TimeDiscretization.GetMinimumBufferSize(bdf), exp_size)

        # test with passing a ProcessInfo
        coeffs = bdf.ComputeBDFCoefficients(model_part.ProcessInfo)
        self.assertEqual(len(coeffs), exp_size)
        for coeff, exp_coeff in zip(coeffs, exp_results):
            self.assertAlmostEqual(coeff, exp_coeff)

    def test_generic_BDF(self):
        delta_time = 0.11
        prev_delta_time = 0.089

        # Check BDF 1
        time_order = 1
        bdf_1 = KM.TimeDiscretization.BDF(time_order)
        exp_results_1 = self._get_expected_results(time_order, delta_time)
        self._test_bdf_time_discretization(bdf_1, exp_results_1, delta_time)

        # Check BDF 2
        time_order = 2
        bdf_2 = KM.TimeDiscretization.BDF(time_order)
        exp_results_2 = self._get_expected_results(time_order, delta_time, prev_delta_time)
        self._test_bdf_time_discretization(bdf_2, exp_results_2, delta_time, prev_delta_time)

        # Check BDF 3 to 6
        for time_order in range(3,7):
            bdf = KM.TimeDiscretization.BDF(time_order)
            exp_results = self._get_expected_results(time_order, delta_time)
            self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF1(self):
        bdf = KM.TimeDiscretization.BDF1()

        delta_time = 0.11
        exp_results = self._get_expected_results(1, delta_time)

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF2(self):
        bdf = KM.TimeDiscretization.BDF2()

        delta_time = 0.11
        prev_delta_time = 0.089
        exp_results = self._get_expected_results(2, delta_time, prev_delta_time)

        self._test_bdf_time_discretization(bdf, exp_results, delta_time, prev_delta_time)

    def test_BDF3(self):
        bdf = KM.TimeDiscretization.BDF3()

        delta_time = 0.11
        exp_results = self._get_expected_results(3, delta_time)

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF4(self):
        bdf = KM.TimeDiscretization.BDF4()

        delta_time = 0.11
        exp_results = self._get_expected_results(4, delta_time)

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF5(self):
        bdf = KM.TimeDiscretization.BDF5()

        delta_time = 0.11
        exp_results = self._get_expected_results(5, delta_time)

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_BDF6(self):
        bdf = KM.TimeDiscretization.BDF6()

        delta_time = 0.11
        exp_results = self._get_expected_results(6, delta_time)

        self._test_bdf_time_discretization(bdf, exp_results, delta_time)

    def test_Newmark(self):
        gen_alpha = KM.TimeDiscretization.Newmark()
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.25)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.5)

        self.assertEqual(KM.TimeDiscretization.GetMinimumBufferSize(gen_alpha), 2)

    def test_Bossak(self):
        gen_alpha = KM.TimeDiscretization.Bossak()
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.4225)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.8)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)

        self.assertEqual(KM.TimeDiscretization.GetMinimumBufferSize(gen_alpha), 2)

    def test_GeneralizedAlpha(self):
        gen_alpha = KM.TimeDiscretization.GeneralizedAlpha()

        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.4225)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.8)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)
        self.assertAlmostEqual(gen_alpha.GetAlphaF(), 0.0)

        self.assertEqual(KM.TimeDiscretization.GetMinimumBufferSize(gen_alpha), 2)

if __name__ == '__main__':
    KratosUnittest.main()
