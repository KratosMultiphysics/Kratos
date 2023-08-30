import KratosMultiphysics

import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestSaving(KratosUnittest.TestCase):

    def test_Saving(self):
        self._execute_saving_test()

if __name__ == '__main__':
    KratosUnittest.main()