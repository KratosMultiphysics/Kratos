# Import python class test
import unittest

# Import python packages
import numpy as np

# Import xmc classes
import xmc.tools as tools


class TestTools(unittest.TestCase):
    def test_normalInverseCDF(self):
        correct_inverse_cdf_values = [
            -2.3263478740408408,
            -1.180559456612439,
            -0.7461851862161866,
            -0.4215776353171568,
            -0.13689839042801627,
            0.13689839042801613,
            0.4215776353171568,
            0.7461851862161862,
            1.180559456612439,
            2.3263478740408408,
        ]
        cdf_values = np.linspace(0.01, 0.99, num=10)
        inverse_cdf_values = [tools.normalInverseCDF(value) for value in cdf_values]
        for i in range(0, len(inverse_cdf_values)):
            self.assertEqual(inverse_cdf_values[i], correct_inverse_cdf_values[i])

    def test_returnInput(self):
        input1 = [[0, 1], [2, 3]]
        input2 = "input2"
        input3 = [3.14]
        output1, output2 = tools.returnInput(input1, input2)
        output3 = tools.returnInput(input3)
        self.assertEqual(output1, input1)
        self.assertEqual(output2, input2)
        self.assertEqual(output3, input3)


if __name__ == "__main__":
    unittest.main()
