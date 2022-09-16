import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestInitialTime(KratosUnittest.TestCase):

    def testIntervalUtilityDefault(self):
        settings = KratosMultiphysics.Parameters()

        # Default interval [0, 1e30]
        small = 1e-12
        interval = KratosMultiphysics.IntervalUtility(settings)
        self.assertTrue(interval.IsInInterval(0.0))
        self.assertFalse(interval.IsInInterval(-small))

        # Custom interval [25, 1e30]
        small = 1e-12
        settings["interval"][0].SetDouble(25.0)
        settings["interval"][1].SetDouble(1000.0)
        interval = KratosMultiphysics.IntervalUtility(settings)
        self.assertTrue(interval.IsInInterval(25.0))
        self.assertFalse(interval.IsInInterval(25.0 - small))

        # Initial interval [0, 0]
        small = 1e-12
        settings["interval"][0].SetDouble(0.0)
        settings["interval"][1].SetDouble(0.0)
        interval = KratosMultiphysics.IntervalUtility(settings)
        self.assertTrue(interval.IsInInterval(0.0))
        self.assertFalse(interval.IsInInterval(small))
        self.assertFalse(interval.IsInInterval(-small))

        # Initial interval [1000, 1000]
        small = 1e-10       # far from 0, bigger threshold
        settings["interval"][0].SetDouble(1000.0)
        settings["interval"][1].SetDouble(1000.0)
        interval = KratosMultiphysics.IntervalUtility(settings)
        self.assertTrue(interval.IsInInterval(1000.0))
        self.assertFalse(interval.IsInInterval(1000.0 - small))
        self.assertFalse(interval.IsInInterval(1000.0 + small))

        # Custom interval [0, 1000]
        small = 1e-10       # bigger interval, higher absolute tolerance
        settings["interval"][0].SetDouble(0.0)
        settings["interval"][1].SetDouble(1000.0)
        interval = KratosMultiphysics.IntervalUtility(settings)
        self.assertTrue(interval.IsInInterval(0.0))
        self.assertTrue(interval.IsInInterval(1.0))
        self.assertTrue(interval.IsInInterval(1000.0))
        self.assertFalse(interval.IsInInterval(-small))
        self.assertFalse(interval.IsInInterval(1000.0 + small))


if __name__ == '__main__':
    KratosUnittest.main()
