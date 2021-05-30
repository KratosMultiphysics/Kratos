import KratosMultiphysics as KM
import numpy as np
from KratosMultiphysics.NeuralNetworkApplication.mask_zeros_process import MaskZerosProcess
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestMaskZeros(KratosUnittest.TestCase):

    def _execute_mask_zeros_test(self):
        # Data inputs
        data_in = np.array([[1, 0, 3, 0, 5],
                   [7, 0, 9, 0, 11],
                   [13, 0, 15, 0, 17]])
        data_out = np.array([[1, 0, 3, 0, 5],
                    [7, 0, 9, 0, 11],
                    [13, 0, 15, 0, 17]])
        # Processed data
        processed_data_in = np.array([[1, 3, 5],
                   [7, 9, 11],
                   [13, 15, 17]])
        processed_data_out = np.array([[1, 3, 5],
                    [7, 9, 11],
                    [13, 15, 17]])

        parameters = KM.Parameters()
        parameters.AddEmptyValue("objective")
        parameters["objective"].SetString("input")
        process = MaskZerosProcess(parameters)
        [data_in, data_out] = process.Preprocess(data_in, data_out)
        self.assertTrue(all(np.ravel(data_in) == np.ravel(processed_data_in)))

        parameters["objective"].SetString("output")
        process = MaskZerosProcess(parameters)
        [data_in, data_out] = process.Preprocess(data_in, data_out)
        self.assertTrue(all(np.ravel(data_out) == np.ravel(processed_data_out)))

    def test_MaskZeros(self):
        self._execute_mask_zeros_test()

if __name__ == '__main__':
    KratosUnittest.main()