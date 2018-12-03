from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTimeDiscretization(KratosUnittest.TestCase):

    def test_BDF1(self):
        bdf = KM.BDF1()

        delta_time = 0.11

        coeffs = bdf.ComputeBDFCoefficients(delta_time)
        self.assertEqual(len(coeffs), 2)
        self.assertAlmostEqual(coeffs[0], 1.0/delta_time)
        self.assertAlmostEqual(coeffs[1], -1.0/delta_time)

        self.assertEqual(KM.GetMinimumBufferSize(bdf), 2)

    def test_BDF2(self):
        bdf = KM.BDF2()

        delta_time = 0.11
        prev_delta_time = 0.089

        coeffs = bdf.ComputeBDFCoefficients(delta_time, prev_delta_time)
        self.assertEqual(len(coeffs), 3)

        rho = prev_delta_time / delta_time;
        time_coeff = 1.0 / (delta_time * rho * rho + delta_time * rho);

        self.assertAlmostEqual(coeffs[0], time_coeff * (rho * rho + 2.0 * rho))
        self.assertAlmostEqual(coeffs[1], -time_coeff * (rho * rho + 2.0 * rho + 1.0))
        self.assertAlmostEqual(coeffs[2], time_coeff)

        self.assertEqual(KM.GetMinimumBufferSize(bdf), 3)

    def test_BDF3(self):
        bdf = KM.BDF3()

        delta_time = 0.11

        coeffs = bdf.ComputeBDFCoefficients(delta_time)
        self.assertEqual(len(coeffs), 4)
        self.assertAlmostEqual(coeffs[0],  11.0/(6.0*delta_time))
        self.assertAlmostEqual(coeffs[1], -18.0/(6.0*delta_time))
        self.assertAlmostEqual(coeffs[2],   9.0/(6.0*delta_time))
        self.assertAlmostEqual(coeffs[3],  -2.0/(6.0*delta_time))

        self.assertEqual(KM.GetMinimumBufferSize(bdf), 4)

    def test_BDF4(self):
        bdf = KM.BDF4()

        delta_time = 0.11

        coeffs = bdf.ComputeBDFCoefficients(delta_time)
        self.assertEqual(len(coeffs), 5)
        self.assertAlmostEqual(coeffs[0],  25.0/(12.0*delta_time))
        self.assertAlmostEqual(coeffs[1], -48.0/(12.0*delta_time))
        self.assertAlmostEqual(coeffs[2],  36.0/(12.0*delta_time))
        self.assertAlmostEqual(coeffs[3], -16.0/(12.0*delta_time))
        self.assertAlmostEqual(coeffs[4],   3.0/(12.0*delta_time))

        self.assertEqual(KM.GetMinimumBufferSize(bdf), 5)

    def test_BDF5(self):
        bdf = KM.BDF5()

        delta_time = 0.11

        coeffs = bdf.ComputeBDFCoefficients(delta_time)
        self.assertEqual(len(coeffs), 6)
        self.assertAlmostEqual(coeffs[0],  137.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[1], -300.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[2],  300.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[3], -200.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[4],   75.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[5],  -12.0/(60.0*delta_time))

        self.assertEqual(KM.GetMinimumBufferSize(bdf), 6)

    def test_BDF6(self):
        bdf = KM.BDF6()

        delta_time = 0.11

        coeffs = bdf.ComputeBDFCoefficients(delta_time)
        self.assertEqual(len(coeffs), 7)
        self.assertAlmostEqual(coeffs[0],  147.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[1], -360.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[2],  450.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[3], -400.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[4],  225.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[5],  -72.0/(60.0*delta_time))
        self.assertAlmostEqual(coeffs[6],   10.0/(60.0*delta_time))

        self.assertEqual(KM.GetMinimumBufferSize(bdf), 7)

    def test_Newmark(self):
        # testing defaults
        gen_alpha = KM.Newmark()
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.25)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.5)

        # testing with input
        newmark_beta = 0.2
        gen_alpha = KM.Newmark(newmark_beta)
        self.assertAlmostEqual(gen_alpha.GetBeta(), newmark_beta)

        # testing with named input
        newmark_beta = 0.15
        gen_alpha = KM.Newmark(newmark_beta=newmark_beta)
        self.assertAlmostEqual(gen_alpha.GetBeta(), newmark_beta)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

    def test_Bossak(self):
        # testing defaults
        gen_alpha = KM.Bossak()
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.2)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.1225)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)

        # testing with input
        newmark_beta = 0.2
        gen_alpha = KM.Bossak(newmark_beta=newmark_beta)
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.2)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.098)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)

        # testing with input
        alpha_m = -0.2
        newmark_beta = 0.2
        gen_alpha = KM.Bossak(alpha_m, newmark_beta)
        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.3)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.128)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), alpha_m)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

    def test_GeneralizedAlpha(self):
        # testing defaults
        gen_alpha = KM.GeneralizedAlpha()

        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.2)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.1225)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), -0.3)
        self.assertAlmostEqual(gen_alpha.GetAlphaF(), 0.0)

        # testing with input
        alpha_m = -0.3
        alpha_f = 0.1
        gen_alpha = KM.GeneralizedAlpha(alpha_m, alpha_f)

        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.1)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.09)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), alpha_m)
        self.assertAlmostEqual(gen_alpha.GetAlphaF(), alpha_f)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

        # testing with input
        alpha_m = -0.3
        alpha_f = 0.1
        newmark_beta = 0.1
        gen_alpha = KM.GeneralizedAlpha(alpha_m, alpha_f, newmark_beta)

        self.assertAlmostEqual(gen_alpha.GetBeta(), 0.1)
        self.assertAlmostEqual(gen_alpha.GetGamma(), 0.036)
        self.assertAlmostEqual(gen_alpha.GetAlphaM(), alpha_m)
        self.assertAlmostEqual(gen_alpha.GetAlphaF(), alpha_f)

        self.assertEqual(KM.GetMinimumBufferSize(gen_alpha), 2)

if __name__ == '__main__':
    KratosUnittest.main()
