# Import python class test
import unittest

# Import python libraries
import json

# Import xmc classes
import xmc.momentEstimator as me
import xmc.methodDefs_momentEstimator.computeCentralMoments as mdccm


class TestMomentEstimator(unittest.TestCase):

    def test_update(self):
        # dimension = 0
        dimension = 0 # len(samples) = 2**dimension
        list_values = [[1.0],[2.0],[3.0],[4.0],[5.0],[6.0]]
        true_power_sums = [[21.0],[91.0],[441.0],[2275.0],[12201.0],[67171.0],[376761.0],[2142595.0],[12313161.0],[71340451.0]]

        # read parameters
        parametersPath = "parameters/parameters_test_momentEstimator.json"
        with open(parametersPath,'r') as parameter_file:
            parameters = json.load(parameter_file)

        for order in [1,2,5]:
            parameters["momentEstimatorInpuctDict"]["order"] = order
            parameters["momentEstimatorInpuctDict"]["updatedPowerSums"] = "xmc.methodDefs_momentEstimator.updatePowerSums.updatePowerSumsOrder"+str(2*order)+"Dimension0" # required order is 2 * order

            # build momentEstimator class
            test_me = me.MomentEstimator(**parameters["momentEstimatorInpuctDict"])

            # update power sums
            for i in range (len(list_values)):
                test_me.update([list_values[i]])

            # test update sample number
            self.assertEqual(test_me._sampleCounter, len(list_values))

            # test update power sums
            for i in range (2*test_me.order):
                self.assertEqual(test_me.powerSums[i], true_power_sums[i])

    def test_value(self):
        # dimension = 0
        dimension = 0 # len(samples) = 2**dimension
        list_values = [[1.0],[2.0],[3.0],[4.0],[5.0]]
        true_estimation_order_1 = 3.0
        true_estimation_order_2 = 2.5
        true_error_order_1 = 0.5
        true_error_order_2 = 0.975

        # read parametrs
        parametersPath = "parameters/parameters_test_momentEstimator.json"
        with open(parametersPath,'r') as parameter_file:
            parameters = json.load(parameter_file)

        # build momentEstimator class
        test_me = me.MomentEstimator(**parameters["momentEstimatorInpuctDict"])

        # update power sums
        test_me.update(list_values)

        # test
        xmc_estimation_order_1 = test_me.value(order=1,isCentral=True,isErrorEstimationRequested=False)
        xmc_estimation_order_2 = test_me.value(order=2,isCentral=True,isErrorEstimationRequested=False)
        xmc_error_order_1 = test_me.value(order=1,isCentral=True,isErrorEstimationRequested=True)
        xmc_error_order_2 = test_me.value(order=2,isCentral=True,isErrorEstimationRequested=True)
        self.assertEqual(xmc_estimation_order_1,true_estimation_order_1)
        self.assertEqual(xmc_estimation_order_2,true_estimation_order_2)
        self.assertEqual(xmc_error_order_1,true_error_order_1)
        self.assertAlmostEqual(xmc_error_order_2,true_error_order_2) # fail if the two objects are unequal
                                                                     # as determined by their difference rounded
                                                                     # to the given number of decimal places(default 7)
                                                                     # (according to unittest docs)

class TestCombinedMomentEstimator(unittest.TestCase):

    def test_update(self):
        # dimension = 0
        dimension = 0 # len(samples) = 2**dimension
        Q=2.0
        list_values = [[[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],[10*Q],[10*Q*Q],10]]

        true_power_sums = [[120.0],[240.0],[120.0],[240.0],[120.0],[240.0],[120.0],[240.0],[120.0],[240.0]]

        # read parameters
        parametersPath = "parameters/parameters_test_combinedMomentEstimator.json"
        with open(parametersPath,'r') as parameter_file:
            parameters = json.load(parameter_file)

        for order in [1,5]:
            parameters["momentEstimatorInpuctDict"]["order"] = order
            # build momentEstimator class
            test_me = me.CombinedMomentEstimator(**parameters["momentEstimatorInpuctDict"])

            # update power sums
            for i in range (len(list_values)):
                test_me.update([[list_values[i]]])

            # test update sample number
            self.assertEqual(test_me._sampleCounter,60)

            # test update power sums
            for i in range (2*test_me.order):
                self.assertEqual(test_me.powerSums[i], true_power_sums[i])

            S1 = test_me.powerSums[0][0]
            S2 = test_me.powerSums[1][0]
            h1 = mdccm.computeCentralMomentsOrderOneDimensionZero(S1,test_me._sampleCounter)
            h2 = mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased(S1,S2,test_me._sampleCounter)

            self.assertEqual(h1,2.0)
            self.assertEqual(h2,0.0)

    def test_value(self):
        # dimension = 0
        dimension = 0 # len(samples) = 2**dimension
        Q=2.0
        list_values = [[[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],10],[[10*Q],[10*Q*Q],10]]

        true_estimation_order_1 = 2.0
        true_estimation_order_2 = 0.0
        true_error_order_1 = 0.0

        # read parameters
        parametersPath = "parameters/parameters_test_combinedMomentEstimator.json"
        with open(parametersPath,'r') as parameter_file:
            parameters = json.load(parameter_file)

        # build momentEstimator class
        test_me = me.CombinedMomentEstimator(**parameters["momentEstimatorInpuctDict"])

        # update power sums
        for i in range (len(list_values)):
            test_me.update([[list_values[i]]])

        # test
        xmc_estimation_order_1 = test_me.value(order=1,isCentral=True,isErrorEstimationRequested=False)
        xmc_estimation_order_2 = test_me.value(order=2,isCentral=True,isErrorEstimationRequested=False)
        xmc_error_order_1 = test_me.value(order=1,isCentral=True,isErrorEstimationRequested=True)
        self.assertEqual(xmc_estimation_order_1,true_estimation_order_1)
        self.assertEqual(xmc_estimation_order_2,true_estimation_order_2)
        self.assertEqual(xmc_error_order_1,true_error_order_1)

if __name__ == '__main__':
    unittest.main()