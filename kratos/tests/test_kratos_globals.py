import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestKratosGlobals(KratosUnittest.TestCase):

    def test_has_variable(self):
        var_name = "VELOCITY"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasVariable(var_name))
        var_name = "VELOCITY2"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasVariable(var_name))

    def test_get_variable_type(self):
        var_name = "TEMPERATURE"
        self.assertEqual(KratosMultiphysics.KratosGlobals.GetVariableType(var_name), "Double")
        var_name = "VELOCITY"
        self.assertEqual(KratosMultiphysics.KratosGlobals.GetVariableType(var_name), "Array")

    def test_get_variable(self):
        var_name = "TEMPERATURE"
        self.assertEqual(KratosMultiphysics.KratosGlobals.GetVariable(var_name), KratosMultiphysics.TEMPERATURE)
        var_name = "VELOCITY"
        self.assertEqual(KratosMultiphysics.KratosGlobals.GetVariable(var_name), KratosMultiphysics.VELOCITY)
        var_name = "VELOCITY2"
        self.assertRaises(Exception, KratosMultiphysics.KratosGlobals.GetVariable, var_name)

    def test_has_flag(self):
        flag_name = "STRUCTURE"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasFlag(flag_name))
        flag_name = "STRUCTURE2"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasFlag(flag_name))

    def test_get_flag(self):
        flag_name = "STRUCTURE"
        KratosMultiphysics.KratosGlobals.GetFlag(flag_name)
        # self.assertEqual(KratosMultiphysics.KratosGlobals.GetFlag(flag_name), KratosMultiphysics.STRUCTURE) # TODO: This should work, they are different objects

    def test_has_constitutive_law(self):
        law_name = "ConstitutiveLaw"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasConstitutiveLaw(law_name))
        law_name = "ConstitutiveLaw2"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasConstitutiveLaw(law_name))

    def test_get_constitutive_law(self):
        law_name = "ConstitutiveLaw"
        KratosMultiphysics.KratosGlobals.GetConstitutiveLaw(law_name)

    def test_has_modeler(self):
        modeler_name = "Modeler"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasModeler(modeler_name))
        modeler_name = "Modeler2"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasModeler(modeler_name))

    def test_get_modeler(self):
        modeler_name = "Modeler"
        KratosMultiphysics.KratosGlobals.GetModeler(modeler_name)

    def test_has_geometry(self):
        geometry_name = "Line2D2"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasGeometry(geometry_name))
        geometry_name = "Lime2D2"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasGeometry(geometry_name))

    def test_get_geometry(self):
        geometry_name = "Line2D2"
        KratosMultiphysics.KratosGlobals.GetGeometry(geometry_name)

    def test_has_condition(self):
        condition_name = "PointCondition3D1N"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasCondition(condition_name))
        condition_name = "PointKondition3D1N"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasCondition(condition_name))

    def test_get_condition(self):
        condition_name = "PointCondition3D1N"
        KratosMultiphysics.KratosGlobals.GetCondition(condition_name)

    def test_has_element(self):
        element_name = "Element3D1N"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasElement(element_name))
        element_name = "Element3D1Ñ"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasElement(element_name))

    def test_get_element(self):
        element_name = "Element3D1N"
        KratosMultiphysics.KratosGlobals.GetElement(element_name)

    def test_has_master_slave_constraint(self):
        constraint_name = "MasterSlaveConstraint"
        self.assertTrue(KratosMultiphysics.KratosGlobals.HasMasterSlaveConstraint(constraint_name))
        constraint_name = "MasterSlaveKonsstraint"
        self.assertFalse(KratosMultiphysics.KratosGlobals.HasMasterSlaveConstraint(constraint_name))

    def test_get_master_slave_constraint(self):
        constraint_name = "MasterSlaveConstraint"
        KratosMultiphysics.KratosGlobals.GetMasterSlaveConstraint(constraint_name)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
