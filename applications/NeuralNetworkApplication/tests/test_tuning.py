import KratosMultiphysics

import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestLoadingHypermodel(KratosUnittest.TestCase):

    def test_LoadingHypermodel(self):
        self._execute_loading_hypermodel_test()

class TestBuildingHypermodel(KratosUnittest.TestCase):

    def test_BuildingHypermodel(self):
        self._execute_building_hypermodel_test()

class TestTuning(KratosUnittest.TestCase):

    def test_Tuning(self):
        self._execute_tuning_test()

if __name__ == '__main__':
    KratosUnittest.main()