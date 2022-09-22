# Core imports
import KratosMultiphysics
from KratosMultiphysics import KratosUnittest

# STD imports
import platform


class TestIntervalUtility(KratosUnittest.TestCase):

    @staticmethod
    def MakeLinspace(begin: float, end: float, number_of_samples: int) -> "list[float]":
        step = (end - begin) / (number_of_samples - 1)
        out = [begin + i_sample * step for i_sample in range(number_of_samples)]
        out[-1] = end
        return out

    @property
    def resolution(self) -> int:
        return 25

    def CheckRange(self,
                   interval: KratosMultiphysics.IntervalUtility,
                   begin: float,
                   end: float,
                   check_boundaries = True,
                   check_below = True,
                   check_above = True) -> None:
        if check_boundaries:
            self.assertAlmostEqual(begin, interval.GetIntervalBegin())
            self.assertAlmostEqual(end, interval.GetIntervalEnd())
        size = end - begin

        # Check range strictly below begin
        if check_below:
            for value in self.MakeLinspace(begin - size, begin, self.resolution)[:-1]:
                self.assertFalse(interval.IsInInterval(value), msg=f"{begin} {value} {end}")

        # Check lower boundary
        self.assertTrue(interval.IsInInterval(begin))

        # Check range strictly in interval
        for value in self.MakeLinspace(begin, end, self.resolution)[1:-1]:
            self.assertTrue(interval.IsInInterval(value))

        # Check upper boundary
        self.assertTrue(interval.IsInInterval(end), msg = f"{end} {interval.GetIntervalEnd()}")

        # Check range strictly above interval
        if check_above:
            for value in self.MakeLinspace(end, end + size, self.resolution)[1:]:
                self.assertFalse(interval.IsInInterval(value), msg=f"{begin} {value} {end}")

    @property
    def min(self) -> float:
        return -1e300

    @property
    def max(self) -> float:
        return 1e300

    def CheckDefaultRange(self, interval: KratosMultiphysics.IntervalUtility) -> None:
        """Same as CheckRange but ranges below and above the interval are not checked (not representable)."""

        # Check lower boundary
        self.assertTrue(interval.IsInInterval(self.min))

        # Check range strictly in interval
        for value in self.MakeLinspace(self.min, self.max, 100 * self.resolution)[1:-1]:
            self.assertTrue(interval.IsInInterval(value))

        # Check upper boundary
        self.assertTrue(interval.IsInInterval(self.max))

    def test_ConstructFromDefaultParameters(self) -> None:
        # Default constructor
        interval = KratosMultiphysics.IntervalUtility()
        self.CheckDefaultRange(interval)

        # Empty Parameters
        parameters = KratosMultiphysics.Parameters()
        interval = KratosMultiphysics.IntervalUtility(parameters)
        self.CheckDefaultRange(interval)

        # Parameters with irrelevant members
        parameters = KratosMultiphysics.Parameters("""{
            "irrelevant" : "setting",
            "useless_array" : [-1,"Begin"]
        }""")
        interval = KratosMultiphysics.IntervalUtility(parameters)
        self.CheckDefaultRange(interval)

    def test_ExactlyRepresentableBoundaries(self) -> None:
        for begin, end in ((1.0, 2.0), (0.5, 4.0), (0.125, 0.25)):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(begin)
            parameters["interval"].Append(end)
            self.CheckRange(KratosMultiphysics.IntervalUtility(parameters), begin, end)

    def test_NotRepresentableBoundaries(self) -> None:
        for begin, end in ((0.3, 2.0), (0.25, 0.3), (1.0/3.0, 0.7)):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(begin)
            parameters["interval"].Append(end)
            self.CheckRange(KratosMultiphysics.IntervalUtility(parameters), begin, end)

    def test_StringBoundaries(self) -> None:
        # Check begin
        for end in (-2.0, -0.3, -1e10):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append("Begin")
            parameters["interval"].Append(end)
            interval = KratosMultiphysics.IntervalUtility(parameters)
            self.CheckRange(interval, self.min, end, check_boundaries=False, check_below=False)

        # Check end
        for begin in (-2.0, -0.3, -1e10):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(begin)
            parameters["interval"].Append("End")
            interval = KratosMultiphysics.IntervalUtility(parameters)
            self.CheckRange(interval, begin, self.max, check_boundaries=False, check_above=False)

        # Check begin & end
        parameters = KratosMultiphysics.Parameters("""{"interval" : ["Begin", "End"]}""")
        self.CheckDefaultRange(KratosMultiphysics.IntervalUtility(parameters))

    def test_DegenerateInterval(self) -> None:
        for value in (-1e10, 0.0, 1e10):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(value)
            parameters["interval"].Append(value)
            self.assertTrue(KratosMultiphysics.IntervalUtility(parameters).IsInInterval(value))

    def test_InvalidParameters(self) -> None:
        invalid_parameters = [KratosMultiphysics.Parameters('{"interval":' + value + '}') for value in [
            "1", # invalid type: numeric
            "{ }", # invalid type: object
            '"[0, 1]"', # invalid type: string
            "[]", # empty array
            "[1]", # incomplete array
            "[1, 2, 3]", # array too long
            '["End", "2"]', # mismplaced "End"
            '[2, "Begin"]', # misplaced "Begin"
            '["You expected Begin", "But it was me, Dio!"]', # invalid strings
            "[[0], [1]]", # invalid types in array: nested array
            '[{"0" : 0}, {"1" : 1}]', # invalid types in array: objects
            "[1, 0]", # invalid order: positives
            "[10, -12]", # invalid order: mixes signs
            "[-2, -5]", # invalid order: negatives
        ]]
        for parameters in invalid_parameters:
            try:
                KratosMultiphysics.IntervalUtility(parameters)
                self.assertTrue(False, msg = f"Expected the following parameters to trigger an exception, but they did not: {parameters}")
            except Exception as exception:
                pass # expected


class TestDiscreteIntervalUtility(KratosUnittest.TestCase):

    @property
    def resolution(self) -> int:
        return 25

    @staticmethod
    def MakeLinspace(begin: int, end: int, number_of_samples: int) -> "list[int]":
        step = (end - begin) / float(number_of_samples - 1)
        out = [round(begin + i_sample * step) for i_sample in range(number_of_samples)]
        out[-1] = end
        return out

    def CheckRange(self,
                   interval: KratosMultiphysics.DiscreteIntervalUtility,
                   begin: int,
                   end: int,
                   check_boundaries = True,
                   check_below = True,
                   check_above = True) -> None:
        if check_boundaries:
            self.assertAlmostEqual(begin, interval.GetIntervalBegin())
            self.assertAlmostEqual(end, interval.GetIntervalEnd())
        size = end - begin

        # Check range strictly below begin
        if check_below:
            for value in self.MakeLinspace(begin - size, begin, self.resolution)[:-1]:
                self.assertFalse(interval.IsInInterval(value), msg=f"{begin} {value} {end}")

        # Check lower boundary
        self.assertTrue(interval.IsInInterval(begin))

        # Check range strictly in interval
        for value in self.MakeLinspace(begin, end, self.resolution)[1:-1]:
            self.assertTrue(interval.IsInInterval(value))

        # Check upper boundary
        self.assertTrue(interval.IsInInterval(end))

        # Check range strictly above interval
        if check_above:
            for value in self.MakeLinspace(end, end + size, self.resolution)[1:]:
                self.assertFalse(interval.IsInInterval(value), msg=f"{begin} {value} {end}")

    @property
    def max(self) -> int:
        """Workaround for python's autogrowing int type."""
        bit_count = int(platform.architecture()[0][:-3])
        unsigned_bit_count = (bit_count - 1) // 2 - 1
        return 1 << unsigned_bit_count

    @property
    def min(self) -> int:
        return -(self.max - 1)

    def CheckDefaultRange(self, interval: KratosMultiphysics.DiscreteIntervalUtility) -> None:
        """Same as CheckRange but ranges below and above the interval are not checked (not representable)."""

        # Check lower boundary
        self.assertTrue(interval.IsInInterval(self.min))

        # Check range strictly in interval
        for value in self.MakeLinspace(self.min, self.max, 100 * self.resolution)[1:-1]:
            self.assertTrue(interval.IsInInterval(value))

        # Check upper boundary
        self.assertTrue(interval.IsInInterval(self.max))

    def test_ConstructFromDefaultParameters(self) -> None:
        # Default constructor
        interval = KratosMultiphysics.DiscreteIntervalUtility()
        self.CheckDefaultRange(interval)

        # Empty Parameters
        parameters = KratosMultiphysics.Parameters()
        interval = KratosMultiphysics.DiscreteIntervalUtility(parameters)
        self.CheckDefaultRange(interval)

        # Parameters with irrelevant members
        parameters = KratosMultiphysics.Parameters("""{
            "irrelevant" : "setting",
            "useless_array" : [-1,"Begin"]
        }""")
        interval = KratosMultiphysics.DiscreteIntervalUtility(parameters)
        self.CheckDefaultRange(interval)

    def test_NumericBoundaries(self) -> None:
        for begin, end in ((1, 30), (-200, 40), (-100, -50)):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(begin)
            parameters["interval"].Append(end)
            self.CheckRange(KratosMultiphysics.DiscreteIntervalUtility(parameters), begin, end)

    def test_StringBoundaries(self) -> None:
        # Begin and end checks cannot be performed correctly in pyhton,
        # because the builtin int type grows its byte representation
        # when approaching its boundaries => it's no longer interpretable
        # as a standard int in C++.
        # Check begin
        for end in (-2, -3, -123456):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append("Begin")
            parameters["interval"].Append(end)
            try:
                interval = KratosMultiphysics.DiscreteIntervalUtility(parameters)
            except Exception as exception:
                print(end, parameters)
                raise exception
            self.CheckRange(interval, self.min, end, check_boundaries=False, check_below=False)

        # Check end
        for begin in (2, 3, 123456):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(begin)
            parameters["interval"].Append("End")
            interval = KratosMultiphysics.DiscreteIntervalUtility(parameters)
            self.CheckRange(interval, begin, self.max, check_boundaries=False, check_above=False)

        # Check begin & end
        parameters = KratosMultiphysics.Parameters("""{"interval" : ["Begin", "End"]}""")
        self.CheckDefaultRange(KratosMultiphysics.DiscreteIntervalUtility(parameters))

    def test_DegenerateInterval(self) -> None:
        for value in (-123456, 0, 123456):
            parameters = KratosMultiphysics.Parameters()
            parameters.AddEmptyArray("interval")
            parameters["interval"].Append(value)
            parameters["interval"].Append(value)
            self.assertTrue(KratosMultiphysics.DiscreteIntervalUtility(parameters).IsInInterval(value))

    def test_InvalidParameters(self) -> None:
        invalid_parameters = [KratosMultiphysics.Parameters('{"interval":' + value + '}') for value in [
            "1", # invalid type: numeric
            "{ }", # invalid type: object
            '"[0, 1]"', # invalid type: string
            "[]", # empty array
            "[1]", # incomplete array
            "[1, 2, 3]", # array too long
            '["End", "2"]', # mismplaced "End"
            '[2, "Begin"]', # misplaced "Begin"
            '["You expected Begin", "But it was me, Dio!"]', # invalid strings
            "[[0], [1]]", # invalid types in array: nested array
            '[{"0" : 0}, {"1" : 1}]', # invalid types in array: objects
            "[1, 0]", # invalid order: positives
            "[10, -12]", # invalid order: mixes signs
            "[-2, -5]", # invalid order: negatives
        ]]
        for parameters in invalid_parameters:
            try:
                KratosMultiphysics.DiscreteIntervalUtility(parameters)
                self.assertTrue(False, msg = f"Expected the following parameters to trigger an exception, but they did not: {parameters}")
            except Exception as exception:
                pass # expected

    def testIntervalUtilityDefault(self):
        settings = KratosMultiphysics.Parameters()

        # Default interval [-inf, inf]
        small = 1e-12
        interval = KratosMultiphysics.IntervalUtility(settings)
        self.assertTrue(interval.IsInInterval(0.0))

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



if __name__ == "__main__":
    KratosUnittest.main()
