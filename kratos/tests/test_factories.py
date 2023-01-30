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

if __name__ == '__main__':
    KratosUnittest.main()
