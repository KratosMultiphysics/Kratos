# Import python class test
import unittest

# Import Python libraries
from json import load
from itertools import product as itprod
import numpy.random as npr

# Import xmc classes
from xmc.randomGeneratorWrapper import NumPyGenerator


class TestNumPyGenerator(unittest.TestCase):
    """Test case for :py:class:`xmc.randomGeneratorWrapper.NumPyGenerator`"""

    def test_init(self):
        """Test constructor method."""
        # Test multiple valid combinations of arguments
        # List valid values
        seeds = (None, 14)
        methods = ("normal", "uniform")
        params = ((), [], [0, 1], (0, 1, 3))
        # Create all combinations of methods and parameters
        distributions = [[(m, p)] for m, p in itprod(methods, params)]
        # Add the list of all distributions (i.e. (m,p) couples above) at the end
        distributions.append([dist[0] for dist in distributions])
        # Test all possible combinations
        for seed, dist in itprod(seeds, distributions):
            with self.subTest(seed=seed, distributions=dist):
                gen = NumPyGenerator(seed=seed, distributions=dist)
        # Test exception for invalid arguments
        distributions = [["normal", 0, 1]]
        with self.subTest(
            msg="Test handling of incompatible 'distributions' arguments",
            distributions=distributions,
        ):
            with self.assertRaises(Exception):
                gen = NumPyGenerator(distributions=distributions)

    def test_realisation(self):
        """Test realisation method."""
        # Test reproducibility
        seed = 2021
        ref_rng = npr.default_rng(seed)
        dist = ("normal", (0, 1, 10))
        gen = NumPyGenerator(seed=seed, distributions=[dist])
        # Generate reference event
        ref_event = getattr(ref_rng, dist[0])(*dist[1])
        ref_event = ref_event.astype(float).tolist()  # convert it
        # Generate event to test
        event = gen.realisation()[0]  # unpack
        # Compare reference and received events
        self.assertEqual(event, ref_event)

    def test_realisation_fromJSON(self):
        "Compare realisations with stored data."
        with open("parameters/parameters_randomGenerator.json", "r") as parameter_file:
            parameters = load(parameter_file)
        expected = parameters["expected"]
        rng = NumPyGenerator(**parameters["randomGeneratorInputDictionary"])
        for i in range(0,len(expected),2):
            self.assertEqual(rng.realisation(), [expected[i],expected[i+1]])


if __name__ == "__main__":
    unittest.main(verbosity=2)
