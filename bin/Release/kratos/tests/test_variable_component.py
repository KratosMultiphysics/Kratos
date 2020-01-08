from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestVariableComponent(KratosUnittest.TestCase):
    def test_GetSourceVariable(self):
        component_variable = KratosMultiphysics.VELOCITY_X
        source_variable = KratosMultiphysics.VELOCITY
        false_source_variable = KratosMultiphysics.DISPLACEMENT
        self.assertEqual(component_variable.GetSourceVariable(), source_variable)
        self.assertNotEqual(component_variable.GetSourceVariable(), false_source_variable)

if __name__ == '__main__':
    KratosUnittest.main()
