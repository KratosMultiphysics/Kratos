import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestFactories(KratosUnittest.TestCase):

    def _auxiliary_test_function_ExplicitBuilder(self, settings, name):
        builder_and_solver = KM.ExplicitBuilderFactory().Create(settings)
        self.assertTrue(KM.ExplicitBuilderFactory().Has(settings["name"].GetString()))
        self.assertEqual(builder_and_solver.Info(), name)

    def test_ExplicitBuilder(self):
        settings = KM.Parameters("""
        {
            "name" : "explicit_builder"
        }
        """)
        self._auxiliary_test_function_ExplicitBuilder(settings, "ExplicitBuilder")

    def _auxiliary_test_function_ExplicitStrategy(self, settings, name):
        this_model = KM.Model()
        model_part = this_model.CreateModelPart("Main")
        strategy = KM.ExplicitStrategyFactory().Create(model_part, settings)
        self.assertTrue(KM.ExplicitStrategyFactory().Has(settings["name"].GetString()))
        self.assertEqual(strategy.Info(), name)

    def test_ExplicitSolvingStrategyRungeKutta4(self):
        settings = KM.Parameters("""
        {
            "name" : "explicit_solving_strategy_runge_kutta_4",
            "explicit_builder_settings" : {
                "name": "explicit_builder"
            }
        }
        """)
        self._auxiliary_test_function_ExplicitStrategy(settings, "ExplicitSolvingStrategyRungeKutta4")

if __name__ == '__main__':
    KratosUnittest.main()
