import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsConditionTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = -86430.42350907
    etalon_value2 = -0.0005
    etalon_value3 = -8643.04235086

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def check_water_pressure(self, test_name, etalon_value):
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(etalon_value, pressure[2])

    def check_displacement(self, test_name, etalon_value):
        file_path = test_helper.get_file_path(os.path.join('test_conditions', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        self.assertAlmostEqual(etalon_value, displacement[2][1])

    def test_PwNormalFluxCondition2D2N(self):
        self.check_water_pressure('test_PwNormalFluxCondition2D2N', self.etalon_value1)

    def test_PwNormalFluxCondition2D3N(self):
        self.check_water_pressure('test_PwNormalFluxCondition2D3N', self.etalon_value1)

    def test_PwNormalFluxCondition2D4N(self):
        self.check_water_pressure('test_PwNormalFluxCondition2D4N', self.etalon_value1)

    def test_PwNormalFluxCondition2D5N(self):
        self.check_water_pressure('test_PwNormalFluxCondition2D5N', self.etalon_value1)

    def test_UPwNormalFluxCondition2D2N(self):
        self.check_water_pressure('test_UPwNormalFluxCondition2D2N', self.etalon_value3)

    def test_UPwNormalFluxCondition2D3N(self):
        self.check_water_pressure('test_UPwNormalFluxCondition2D3N', self.etalon_value3)

    def test_UPwNormalFluxCondition2D4N(self):
        self.check_water_pressure('test_UPwNormalFluxCondition2D4N', self.etalon_value3)

    def test_UPwNormalFluxCondition2D5N(self):
        self.check_water_pressure('test_UPwNormalFluxCondition2D5N', self.etalon_value3)

    def test_UPwFaceLoadCondition2D2N(self):
        self.check_displacement('test_UPwFaceLoadCondition2D2N', self.etalon_value2)

    def test_UPwFaceLoadCondition2D3N(self):
        self.check_displacement('test_UPwFaceLoadCondition2D3N', self.etalon_value2)

    def test_UPwFaceLoadCondition2D4N(self):
        self.check_displacement('test_UPwFaceLoadCondition2D4N', self.etalon_value2)

    def test_UPwFaceLoadCondition2D5N(self):
        self.check_displacement('test_UPwFaceLoadCondition2D5N', self.etalon_value2)

    def test_LineLoadDiffOrderCondition2D3N(self):
        self.check_displacement('test_LineLoadDiffOrderCondition2D3N', self.etalon_value2)

    def test_LineLoadDiffOrderCondition2D4N(self):
        self.check_displacement('test_LineLoadDiffOrderCondition2D4N', self.etalon_value2)

    def test_LineLoadDiffOrderCondition2D5N(self):
        self.check_displacement('test_LineLoadDiffOrderCondition2D5N', self.etalon_value2)

    def test_UPwNormalFaceLoadCondition2D2N(self):
        self.check_displacement('test_UPwNormalFaceLoadCondition2D2N', self.etalon_value2)

    def test_UPwNormalFaceLoadCondition2D3N(self):
        self.check_displacement('test_UPwNormalFaceLoadCondition2D3N', self.etalon_value2)

    def test_UPwNormalFaceLoadCondition2D4N(self):
        self.check_displacement('test_UPwNormalFaceLoadCondition2D4N', self.etalon_value2)

    def test_UPwNormalFaceLoadCondition2D5N(self):
        self.check_displacement('test_UPwNormalFaceLoadCondition2D5N', self.etalon_value2)

    def test_LineNormalLoadDiffOrderCondition2D3N(self):
        self.check_displacement('test_LineNormalLoadDiffOrderCondition2D3N', self.etalon_value2)

    def test_LineNormalLoadDiffOrderCondition2D4N(self):
        self.check_displacement('test_LineNormalLoadDiffOrderCondition2D4N', self.etalon_value2)

    def test_LineNormalLoadDiffOrderCondition2D5N(self):
        self.check_displacement('test_LineNormalLoadDiffOrderCondition2D5N', self.etalon_value2)

    def test_AxisymmetricUPwNormalFaceLoadCondition2D2N(self):
        self.check_displacement('test_AxisymmetricUPwNormalFaceLoadCondition2D2N', self.etalon_value2)

    def test_AxisymmetricUPwNormalFaceLoadCondition2D3N(self):
        self.check_displacement('test_AxisymmetricUPwNormalFaceLoadCondition2D3N', self.etalon_value2)

    def test_AxisymmetricUPwNormalFaceLoadCondition2D4N(self):
        self.check_displacement('test_AxisymmetricUPwNormalFaceLoadCondition2D4N', self.etalon_value2)

    def test_AxisymmetricUPwNormalFaceLoadCondition2D5N(self):
        self.check_displacement('test_AxisymmetricUPwNormalFaceLoadCondition2D5N', self.etalon_value2)
        
if __name__ == '__main__':
    KratosUnittest.main()
