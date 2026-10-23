import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestEntityDataContainer(KratosUnittest.TestCase):
    """Tests for the Kratos::Future entity data containers (DataContainer-backed parallel storage)."""

    def _create_model_part(self, model, buffer_size=1):
        model_part = model.CreateModelPart("test")
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        model_part.SetBufferSize(buffer_size)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        return model_part

    def test_registration(self):
        edc = KM.Future.EntityDataContainer(1, 4)
        self.assertEqual(edc.NumberOfEntities(), 0)
        self.assertEqual(edc.RegisterEntityId(11), 0)
        self.assertEqual(edc.RegisterEntityId(7), 1)
        self.assertEqual(edc.RegisterEntityId(11), 0)  # idempotent
        self.assertTrue(edc.HasEntity(7))
        self.assertEqual(edc.Index(7), 1)
        self.assertEqual(edc.NumberOfEntities(), 2)

        with self.assertRaises(RuntimeError):
            edc.Index(99)

    def test_values_and_growth(self):
        edc = KM.Future.EntityDataContainer(1, 2)  # tiny chunk size to force growth
        edc.AddVariable(KM.PRESSURE)

        for entity_id in range(1, 6):
            edc.RegisterEntityId(entity_id)
            edc.SetValue(entity_id, KM.PRESSURE, 1.5 * entity_id)

        # Growth preserved the previously written values
        for entity_id in range(1, 6):
            self.assertEqual(edc.GetValue(entity_id, KM.PRESSURE), 1.5 * entity_id)

    def test_numpy_span_access(self):
        edc = KM.Future.EntityDataContainer(1, 4)
        for entity_id in [10, 20, 30]:
            edc.RegisterEntityId(entity_id)
        accessor = edc.AddVariable(KM.PRESSURE)

        span = edc.GetDataContainer().GetDataSpan(accessor)
        self.assertEqual(len(span), edc.Capacity())
        span[edc.Index(20)] = 4.5

        self.assertEqual(edc.GetValue(20, KM.PRESSURE), 4.5)

    def test_model_part_container_nodes(self):
        model = KM.Model()
        model_part = self._create_model_part(model, buffer_size=3)

        mpdc = KM.Future.ModelPartDataContainer(model_part)
        nodes_data = mpdc.Nodes()
        self.assertEqual(nodes_data.NumberOfEntities(), 3)
        self.assertEqual(nodes_data.GetBufferSize(), 3)

        nodes_data.AddHistoricalVariable(KM.TEMPERATURE)
        nodes_data.AddVariable(KM.PRESSURE)

        for node in model_part.Nodes:
            nodes_data.SetValue(node.Id, KM.TEMPERATURE, 100.0 + node.Id)

        # Advance the buffer in lockstep with the legacy workflow
        model_part.CloneTimeStep(1.0)
        mpdc.CloneStepData(KM.StepCategory.TimeStep)
        nodes_data.SetValue(1, KM.TEMPERATURE, 111.0)

        self.assertEqual(nodes_data.GetValue(1, KM.TEMPERATURE), 111.0)
        self.assertEqual(nodes_data.GetValue(1, KM.TEMPERATURE, 1), 101.0)  # previous step
        self.assertEqual(nodes_data.GetValue(2, KM.TEMPERATURE), 102.0)     # cloned

    def test_no_cross_contamination(self):
        model = KM.Model()
        model_part = self._create_model_part(model, buffer_size=2)
        node = model_part.GetNode(1)

        # Seed the LEGACY storage first
        node.SetSolutionStepValue(KM.TEMPERATURE, 25.0)
        node.SetValue(KM.PRESSURE, 5.0)

        mpdc = KM.Future.ModelPartDataContainer(model_part)
        nodes_data = mpdc.Nodes()
        nodes_data.AddHistoricalVariable(KM.TEMPERATURE)
        nodes_data.AddVariable(KM.PRESSURE)

        # New path starts at zero; writes do not leak in either direction
        self.assertEqual(nodes_data.GetValue(1, KM.TEMPERATURE), 0.0)
        nodes_data.SetValue(1, KM.TEMPERATURE, 999.0)
        nodes_data.SetValue(1, KM.PRESSURE, 888.0)
        self.assertEqual(node.GetSolutionStepValue(KM.TEMPERATURE), 25.0)
        self.assertEqual(node.GetValue(KM.PRESSURE), 5.0)

        node.SetSolutionStepValue(KM.TEMPERATURE, 26.0)
        node.SetValue(KM.PRESSURE, 6.0)
        self.assertEqual(nodes_data.GetValue(1, KM.TEMPERATURE), 999.0)
        self.assertEqual(nodes_data.GetValue(1, KM.PRESSURE), 888.0)

    def test_update_after_new_entities(self):
        model = KM.Model()
        model_part = self._create_model_part(model)

        mpdc = KM.Future.ModelPartDataContainer(model_part)
        mpdc.Nodes().AddVariable(KM.PRESSURE)
        mpdc.Nodes().SetValue(1, KM.PRESSURE, 1.5)

        model_part.CreateNewNode(4, 1.0, 1.0, 0.0)
        self.assertFalse(mpdc.Nodes().HasEntity(4))

        mpdc.Update()

        self.assertTrue(mpdc.Nodes().HasEntity(4))
        self.assertEqual(mpdc.Nodes().GetValue(4, KM.PRESSURE), 0.0)
        self.assertEqual(mpdc.Nodes().GetValue(1, KM.PRESSURE), 1.5)  # preserved across growth

    def test_elements_conditions_geometries(self):
        model = KM.Model()
        model_part = model.CreateModelPart("test")
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        properties = model_part.GetProperties()[1]
        elem1 = model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
        elem2 = model_part.CreateNewElement("Element2D3N", 2, [1, 2, 3], properties)
        cond = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewGeometry("Triangle2D3", 1, [1, 2, 3])

        mpdc = KM.Future.ModelPartDataContainer(model_part)
        self.assertEqual(mpdc.Elements().NumberOfEntities(), 2)
        self.assertEqual(mpdc.Conditions().NumberOfEntities(), 1)
        self.assertEqual(mpdc.Geometries().NumberOfEntities(), 1)

        # Element values are independent per Id through the new path
        elements_data = mpdc.Elements()
        elements_data.AddVariable(KM.PRESSURE)
        elements_data.SetValue(elem1.Id, KM.PRESSURE, 1.0)
        elements_data.SetValue(elem2.Id, KM.PRESSURE, 2.0)
        self.assertEqual(elements_data.GetValue(elem1.Id, KM.PRESSURE), 1.0)
        self.assertEqual(elements_data.GetValue(elem2.Id, KM.PRESSURE), 2.0)

        # No cross-contamination with the legacy element storage
        elem1.SetValue(KM.PRESSURE, 10.0)
        self.assertEqual(elem1.GetValue(KM.PRESSURE), 10.0)
        self.assertEqual(elements_data.GetValue(elem1.Id, KM.PRESSURE), 1.0)

        conditions_data = mpdc.Conditions()
        conditions_data.AddVariable(KM.TEMPERATURE)
        conditions_data.SetValue(cond.Id, KM.TEMPERATURE, 42.0)
        self.assertEqual(conditions_data.GetValue(cond.Id, KM.TEMPERATURE), 42.0)
        self.assertFalse(cond.Has(KM.TEMPERATURE))  # legacy storage untouched

        geometries_data = mpdc.Geometries()
        geometries_data.AddVariable(KM.DISTANCE)
        geometries_data.SetValue(1, KM.DISTANCE, 3.5)
        self.assertEqual(geometries_data.GetValue(1, KM.DISTANCE), 3.5)

    def test_master_slave_constraints(self):
        model = KM.Model()
        model_part = model.CreateModelPart("test")
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        n1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        KM.VariableUtils().AddDof(KM.PRESSURE, model_part)
        constraint = model_part.CreateNewMasterSlaveConstraint(
            "LinearMasterSlaveConstraint", 1, n1, KM.PRESSURE, n2, KM.PRESSURE, 0.5, 0.0)

        mpdc = KM.Future.ModelPartDataContainer(model_part)
        constraints_data = mpdc.MasterSlaveConstraints()
        self.assertEqual(constraints_data.NumberOfEntities(), 1)

        constraints_data.AddVariable(KM.DISTANCE)
        constraints_data.SetValue(constraint.Id, KM.DISTANCE, 12.5)
        self.assertEqual(constraints_data.GetValue(constraint.Id, KM.DISTANCE), 12.5)

        # No cross-contamination with the constraint's legacy storage
        self.assertFalse(constraint.Has(KM.DISTANCE))
        constraint.SetValue(KM.DISTANCE, 99.0)
        self.assertEqual(constraint.GetValue(KM.DISTANCE), 99.0)
        self.assertEqual(constraints_data.GetValue(constraint.Id, KM.DISTANCE), 12.5)

    def test_errors(self):
        edc = KM.Future.EntityDataContainer(1, 4)
        edc.RegisterEntityId(1)
        edc.AddVariable(KM.PRESSURE)

        with self.assertRaises(RuntimeError):
            edc.GetValue(1, KM.TEMPERATURE)         # un-added variable
        with self.assertRaises(RuntimeError):
            edc.GetValue(9, KM.PRESSURE)            # unknown entity
        with self.assertRaises(RuntimeError):
            edc.GetValue(1, KM.PRESSURE, 1)         # step request on non-historical
        self.assertTrue(edc.Has(KM.PRESSURE))
        self.assertFalse(edc.Has(KM.TEMPERATURE))


if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
