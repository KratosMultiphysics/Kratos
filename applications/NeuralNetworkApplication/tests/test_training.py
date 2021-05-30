import KratosMultiphysics

import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestTraining(KratosUnittest.TestCase):

    def test_Training(self):
        self._execute_training_test()

if __name__ == '__main__':
    KratosUnittest.main()