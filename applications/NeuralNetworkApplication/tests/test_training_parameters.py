import KratosMultiphysics

import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestOptimizer(KratosUnittest.TestCase):

    def test_KerasOptimizer(self):
        self._execute_optimizer_test(optimizer = 'adam')

class TestLossFunction(KratosUnittest.TestCase):

    def test_KerasLoss(self):
        self._execute_loss_function_test(loss = 'rmse')

class TestMetrics(KratosUnittest.TestCase):

    def test_KerasMetric(self):
        self._execute_metric_test(metric = 'rmse', module = 'keras')

class TestCallbacks(KratosUnittest.TestCase):

    def test_EarlyStopping(self):
        self._execute_callback_test(callback = 'early_stopping')

if __name__ == '__main__':
    KratosUnittest.main()