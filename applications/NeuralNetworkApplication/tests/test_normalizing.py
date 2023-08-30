import KratosMultiphysics as KM
import numpy as np
from KratosMultiphysics.NeuralNetworkApplication.normalization_process import NormalizationProcess
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestNormalizing(KratosUnittest.TestCase):
    def _execute_normalizing_test(self):    
        # Data inputs
        data_in = np.array([[0, 0],
                            [0, 4],
                            [4, 0],
                            [4, 4]])
        data_out = np.array([[0, 0],
                            [0, 4],
                            [4, 0],
                            [4, 4]])
        # Processed data
        processed_data_in = np.array([[-1, -1],
                            [-1, 1],
                            [1, -1],
                            [1, 1]])
        processed_data_out = np.array([[-1, -1],
                            [-1, 1],
                            [1, -1],
                            [1, 1]])

        parameters = KM.Parameters()
        parameters.AddEmptyValue("objective")
        parameters["objective"].SetString("input")
        process = NormalizationProcess(parameters)
        [data_in, data_out] = process.Preprocess(data_in, data_out)
        self.assertTrue(all(np.ravel(data_in) == np.ravel(processed_data_in)))

        parameters["objective"].SetString("output")
        process = NormalizationProcess(parameters)
        [data_in, data_out] = process.Preprocess(data_in, data_out)
        self.assertTrue(all(np.ravel(data_out) == np.ravel(processed_data_out)))
    def test_Normalizing(self):
        self._execute_normalizing_test()

if __name__ == '__main__':
    KratosUnittest.main()