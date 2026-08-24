import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.impose_nodal_young_modulus_process import (
    Factory,
)


class TestImposeNodalYoungModulusProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.previous_threads = KratosMultiphysics.ParallelUtilities.GetNumThreads()

    def tearDown(self):
        KratosMultiphysics.ParallelUtilities.SetNumThreads(self.previous_threads)

    @staticmethod
    def _CreateModelPart(num_nodes):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_YOUNG_MODULUS)
        for i in range(num_nodes):
            model_part.CreateNewNode(i + 1, float(i), float(i % 3), 0.25 * i)
        return model, model_part

    @staticmethod
    def _CreateSettings():
        return KratosMultiphysics.Parameters(
            """
            {
                "Parameters" : {
                    "model_part_name" : "main",
                    "variable_name"   : "NODAL_YOUNG_MODULUS",
                    "Young_Modulus_1" : 1.0,
                    "Young_Modulus_2" : 2.0,
                    "Young_Modulus_3" : 3.0,
                    "Young_Modulus_4" : 4.0
                }
            }
            """
        )

    @staticmethod
    def _ExpectedValue(node):
        return 1.0 + 2.0 * node.X + 3.0 * node.Y + 4.0 * node.Z

    @staticmethod
    def _NodalValues(model_part):
        return [
            node.GetSolutionStepValue(KratosDam.NODAL_YOUNG_MODULUS)
            for node in model_part.Nodes
        ]

    def test_factory_returns_valid_process(self):
        model, _ = self._CreateModelPart(4)
        process = Factory(self._CreateSettings(), model)
        self.assertTrue(isinstance(process, KratosMultiphysics.Process))

    def test_execution_assigns_values_and_does_not_create_dof(self):
        model, model_part = self._CreateModelPart(8)
        process = Factory(self._CreateSettings(), model)

        process.ExecuteBeforeSolutionLoop()

        for node in model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KratosDam.NODAL_YOUNG_MODULUS),
                self._ExpectedValue(node),
            )
            # NODAL_YOUNG_MODULUS is a prescribed nodal material field, not a
            # degree of freedom; no DOF may be created by the process.
            self.assertFalse(node.HasDofFor(KratosDam.NODAL_YOUNG_MODULUS))

    def test_multithreaded_execution_matches_single_thread(self):
        # Single-threaded reference
        KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
        model_ref, model_part_ref = self._CreateModelPart(64)
        Factory(self._CreateSettings(), model_ref).ExecuteBeforeSolutionLoop()
        reference_values = self._NodalValues(model_part_ref)

        # Multi-threaded run of the same production callback
        for threads in (2, 4):
            KratosMultiphysics.ParallelUtilities.SetNumThreads(threads)
            model, model_part = self._CreateModelPart(64)
            Factory(self._CreateSettings(), model).ExecuteBeforeSolutionLoop()
            values = self._NodalValues(model_part)
            self.assertEqual(len(values), len(reference_values))
            for value, reference_value in zip(values, reference_values):
                self.assertAlmostEqual(value, reference_value)

    def test_initialize_solution_step_multithreaded(self):
        KratosMultiphysics.ParallelUtilities.SetNumThreads(4)
        model, model_part = self._CreateModelPart(64)
        process = Factory(self._CreateSettings(), model)

        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KratosDam.NODAL_YOUNG_MODULUS),
                self._ExpectedValue(node),
            )
            self.assertFalse(node.HasDofFor(KratosDam.NODAL_YOUNG_MODULUS))

    def _NodalYoungTetraModel(self, nodal_young):
        """Build a single-tetra model with the given nodal NODAL_YOUNG_MODULUS."""
        model, model_part = self._CreateModelPart(4)
        nodes = list(model_part.Nodes)
        coordinates = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        for node, coords in zip(nodes, coordinates):
            node.X = coords[0]
            node.Y = coords[1]
            node.Z = coords[2]
        for node, young in zip(nodes, nodal_young):
            node.SetSolutionStepValue(KratosDam.NODAL_YOUNG_MODULUS, 0, young)
        properties = model_part.GetProperties()[1]
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        law = KratosDam.LinearElastic3DLawNodal()
        geometry = KratosMultiphysics.Tetrahedra3D4(*nodes)
        return model_part, properties, law, geometry

    def _ConstitutiveResponse(self, model_part, properties, law, geometry):
        """Evaluate the law at the tetra centroid for a fixed strain state."""
        shape = KratosMultiphysics.Vector([0.25, 0.25, 0.25, 0.25])
        law.InitializeMaterial(properties, geometry, shape)

        strain_size = law.GetStrainSize()
        strain_vector = KratosMultiphysics.Vector(strain_size)
        stress_vector = KratosMultiphysics.Vector(strain_size)
        constitutive_matrix = KratosMultiphysics.Matrix(strain_size, strain_size)
        for i in range(strain_size):
            strain_vector[i] = 0.0
            stress_vector[i] = 0.0
        strain_vector[0] = 1.0e-3  # uniaxial strain along x

        options = KratosMultiphysics.Flags()
        options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)

        params = KratosMultiphysics.ConstitutiveLawParameters()
        params.SetOptions(options)
        params.SetDeterminantF(1.0)
        params.SetStrainVector(strain_vector)
        params.SetStressVector(stress_vector)
        params.SetConstitutiveMatrix(constitutive_matrix)
        params.SetShapeFunctionsValues(shape)
        params.SetProcessInfo(model_part.ProcessInfo)
        params.SetMaterialProperties(properties)
        params.SetElementGeometry(geometry)
        law.CalculateMaterialResponseCauchy(params)
        return params.GetStressVector()[0], params.GetConstitutiveMatrix()[0, 0]

    def test_constitutive_response_uses_interpolated_young_modulus(self):
        # The constitutive response must be linear in the interpolated nodal
        # Young modulus: rescaling all nodal values by a factor rescales the
        # response by the same factor (the nodal field drives the material).
        model_part_a, props_a, law_a, geom_a = self._NodalYoungTetraModel(
            [10.0, 20.0, 30.0, 40.0]
        )
        model_part_b, props_b, law_b, geom_b = self._NodalYoungTetraModel(
            [20.0, 40.0, 60.0, 80.0]
        )

        stress_a, c00_a = self._ConstitutiveResponse(model_part_a, props_a, law_a, geom_a)
        stress_b, c00_b = self._ConstitutiveResponse(model_part_b, props_b, law_b, geom_b)

        # Doubling every nodal value doubles the interpolated E, hence doubles
        # the linear-elastic stress and constitutive matrix.
        self.assertAlmostEqual(stress_b, 2.0 * stress_a, 12)
        self.assertAlmostEqual(c00_b, 2.0 * c00_a, 12)

        # The response is deterministic (repeated evaluation is identical).
        stress_a2, c00_a2 = self._ConstitutiveResponse(model_part_a, props_a, law_a, geom_a)
        self.assertAlmostEqual(stress_a2, stress_a, 12)
        self.assertAlmostEqual(c00_a2, c00_a, 12)

    def test_constitutive_response_matches_nodal_interpolation_change(self):
        # A non-uniform nodal field produces the expected interpolated value:
        # changing only some nodal values must change the response, i.e. the
        # response follows the nodal interpolation, not a single material
        # constant.
        model_part_a, props_a, law_a, geom_a = self._NodalYoungTetraModel(
            [10.0, 20.0, 30.0, 40.0]
        )
        model_part_b, props_b, law_b, geom_b = self._NodalYoungTetraModel(
            [10.0, 20.0, 30.0, 80.0]  # only the 4th node doubled
        )

        stress_a, c00_a = self._ConstitutiveResponse(model_part_a, props_a, law_a, geom_a)
        stress_b, c00_b = self._ConstitutiveResponse(model_part_b, props_b, law_b, geom_b)

        # Changing a single nodal value changes the interpolated response.
        self.assertNotEqual(stress_b, stress_a)
        self.assertNotEqual(c00_b, c00_a)

        # With N_i = 1/4, the interpolated E increases from (10+20+30+40)/4=25
        # to (10+20+30+80)/4=35, a factor 35/25 = 1.4 for the elastic response.
        self.assertAlmostEqual(c00_b, (35.0 / 25.0) * c00_a, 12)
        self.assertAlmostEqual(stress_b, (35.0 / 25.0) * stress_a, 12)


if __name__ == "__main__":
    KratosUnittest.main()
