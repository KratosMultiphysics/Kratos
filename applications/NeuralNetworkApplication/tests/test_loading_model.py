import KratosMultiphysics

import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestLoadingModel(KratosUnittest.TestCase):




    def test_LoadingModel(self):
        self._execute_loading_model_test(model = 'model')

if __name__ == '__main__':
    KratosUnittest.main()