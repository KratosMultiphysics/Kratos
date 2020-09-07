import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestFactories(KratosUnittest.TestCase):

    def _auxiliary_test_function_BuilderAndSolver(self, settings, name):
        linear_solver = None
        builder_and_solver = KM.BuilderAndSolverFactory().Create(linear_solver, settings)
        self.assertTrue(KM.BuilderAndSolverFactory().Has(settings["name"].GetString()))
        self.assertEqual(builder_and_solver.Info(), name)

    def test_ResidualBasedEliminationBuilderAndSolver(self):
        settings = KM.Parameters("""
        {
            "name" : "elimination_builder_and_solver"
        }
        """)
        self._auxiliary_test_function_BuilderAndSolver(settings, "ResidualBasedEliminationBuilderAndSolver")

    def test_ResidualBasedEliminationBuilderAndSolverWithConstraints(self):
        settings = KM.Parameters("""
        {
            "name" : "elimination_builder_and_solver_with_constraints"
        }
        """)
        self._auxiliary_test_function_BuilderAndSolver(settings, "ResidualBasedEliminationBuilderAndSolverWithConstraints")

    def test_ResidualBasedBlockBuilderAndSolver(self):
        settings = KM.Parameters("""
        {
            "name" : "block_builder_and_solver"
        }
        """)
        self._auxiliary_test_function_BuilderAndSolver(settings, "ResidualBasedBlockBuilderAndSolver")

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

    def _auxiliary_test_function_ConvergenceCriteria(self, settings, name):
        conv_crit = KM.ConvergenceCriteriaFactory().Create(settings)
        self.assertTrue(KM.ConvergenceCriteriaFactory().Has(settings["name"].GetString()))
        self.assertEqual(conv_crit.Info(), name)

    def test_DisplacementCriteria(self):
        settings = KM.Parameters("""
        {
            "name" : "displacement_criteria"
        }
        """)
        self._auxiliary_test_function_ConvergenceCriteria(settings, "DisplacementCriteria")

    def test_ResidualCriteria(self):
        settings = KM.Parameters("""
        {
            "name" : "residual_criteria"
        }
        """)
        self._auxiliary_test_function_ConvergenceCriteria(settings, "ResidualCriteria")

    def test_And_Criteria(self):
        settings = KM.Parameters("""
        {
            "name" : "and_criteria"
        }
        """)
        self._auxiliary_test_function_ConvergenceCriteria(settings, "And_Criteria")

    def test_Or_Criteria(self):
        settings = KM.Parameters("""
        {
            "name" : "or_criteria"
        }
        """)
        self._auxiliary_test_function_ConvergenceCriteria(settings, "Or_Criteria")

    def test_MixedGenericCriteria(self):
        settings = KM.Parameters("""
        {
            "name" : "mixed_generic_criteria"
        }
        """)
        self._auxiliary_test_function_ConvergenceCriteria(settings, "MixedGenericCriteria")

    def _auxiliary_test_function_Scheme(self, settings, name):
        scheme = KM.SchemeFactory().Create(settings)
        self.assertTrue(KM.SchemeFactory().Has(settings["name"].GetString()))
        self.assertEqual(scheme.Info(), name)

    def test_ResidualBasedIncrementalUpdateStaticScheme(self):
        settings = KM.Parameters("""
        {
            "name" : "static_scheme"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedIncrementalUpdateStaticScheme")

    def test_ResidualBasedIncrementalUpdateStaticSchemeSlip(self):
        settings = KM.Parameters("""
        {
            "name" : "static_slip_scheme"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedIncrementalUpdateStaticSchemeSlip")

    def test_ResidualBasedBossakDisplacementScheme(self):
        settings = KM.Parameters("""
        {
            "name" : "bossak_scheme"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedBossakDisplacementScheme")

    def test_ResidualBasedNewmarkDisplacementScheme(self):
        settings = KM.Parameters("""
        {
            "name" : "newmark_scheme"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedNewmarkDisplacementScheme")

    def test_ResidualBasedPseudoStaticDisplacementScheme(self):
        settings = KM.Parameters("""
        {
            "name" : "pseudo_static_scheme", "rayleigh_beta_variable" : "PRESSURE"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedPseudoStaticDisplacementScheme")

    def test_ResidualBasedBDFDisplacementScheme(self):
        settings = KM.Parameters("""
        {
            "name" : "bdf_displacement_scheme"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedBDFDisplacementScheme")

    def test_ResidualBasedBDFCustomScheme(self):
        settings = KM.Parameters("""
        {
            "name" : "bdf_scheme"
        }
        """)
        self._auxiliary_test_function_Scheme(settings, "ResidualBasedBDFCustomScheme")

    def _auxiliary_test_function_Strategy(self, settings, name):
        this_model = KM.Model()
        model_part = this_model.CreateModelPart("Main")
        strategy = KM.StrategyFactory().Create(model_part, settings)
        self.assertTrue(KM.StrategyFactory().Has(settings["name"].GetString()))
        self.assertEqual(strategy.Info(), name)

    def test_ResidualBasedLinearStrategy(self):
        settings = KM.Parameters("""
        {
            "name" : "linear_strategy",
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "scheme_settings" : {
                "name" : "static_scheme"
            },
            "builder_and_solver_settings" : {
                "name" : "elimination_builder_and_solver"
            }
        }
        """)
        self._auxiliary_test_function_Strategy(settings, "ResidualBasedLinearStrategy")

    def test_ResidualBasedNewtonRaphsonStrategy(self):
        settings = KM.Parameters("""
        {
            "name" : "newton_raphson_strategy",
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "scheme_settings" : {
                "name" : "static_scheme"
            },
            "convergence_criteria_settings" : {
                "name" : "displacement_criteria"
            },
            "builder_and_solver_settings" : {
                "name" : "elimination_builder_and_solver"
            }
        }
        """)
        self._auxiliary_test_function_Strategy(settings, "ResidualBasedNewtonRaphsonStrategy")

    def test_AdaptiveResidualBasedNewtonRaphsonStrategy(self):
        settings = KM.Parameters("""
        {
            "name" : "adaptative_newton_raphson_strategy",
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "scheme_settings" : {
                "name" : "static_scheme"
            },
            "convergence_criteria_settings" : {
                "name" : "displacement_criteria"
            },
            "builder_and_solver_settings" : {
                "name" : "elimination_builder_and_solver"
            }
        }
        """)
        self._auxiliary_test_function_Strategy(settings, "AdaptiveResidualBasedNewtonRaphsonStrategy")

    def test_LineSearchStrategy(self):
        settings = KM.Parameters("""
        {
            "name" : "line_search_strategy",
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "scheme_settings" : {
                "name" : "static_scheme"
            },
            "convergence_criteria_settings" : {
                "name" : "displacement_criteria"
            },
            "builder_and_solver_settings" : {
                "name" : "elimination_builder_and_solver"
            }
        }
        """)
        self._auxiliary_test_function_Strategy(settings, "LineSearchStrategy")

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
